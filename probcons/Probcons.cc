/////////////////////////////////////////////////////////////////
// Main.cc
//
// Main routines for PROBCONS program.
/////////////////////////////////////////////////////////////////

#include "Python.h"


typedef struct
{
    PyObject_HEAD
} Probcons;

Probcons* newProbcons(void);

#include "SafeVector.h"
#include "MultiSequence.h"
#include "Defaults.h"
#include "ScoreType.h"
#include "ProbabilisticModel.h"
#include "EvolutionaryTree.h"
#include "SparseMatrix.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>


static PyObject *ProbconsError;

string parametersInputFilename = "";
string parametersOutputFilename = "no training";
string annotationFilename = "pyconf";

bool enableTraining = false;
bool enableVerbose = false;
bool enableAllPairs = false;
bool enableAnnotation = false;
bool enableViterbi = false;
bool enableClustalWOutput = false;
bool enableTrainEmissions = false;
bool enableAlignOrder = false;

int numConsistencyReps = 2;
int consgap = 0; //whether the consistency transformation should consider gaps in Z, default off=0 

int numPreTrainingReps = 0;
int numIterativeRefinementReps = 100;

float cutoff = 0;
float gapOpenPenalty = 0;
float gapContinuePenalty = 0;

VF initDistrib (NumMatrixTypes);
VF gapOpen (2*NumInsertStates);
VF gapExtend (2*NumInsertStates);
VVF emitPairs (256, VF (256, 1e-10));
VF emitSingle (256, 1e-5);

string alphabet = alphabetDefault;

const int MIN_PRETRAINING_REPS = 0;
const int MAX_PRETRAINING_REPS = 20;
const int MIN_CONSISTENCY_REPS = 0;
const int MAX_CONSISTENCY_REPS = 5;
const int MIN_ITERATIVE_REFINEMENT_REPS = 0;
const int MAX_ITERATIVE_REFINEMENT_REPS = 1000;

/////////////////////////////////////////////////////////////////
// Function prototypes
/////////////////////////////////////////////////////////////////

void PrintHeading();
void PrintParameters (const char *message, const VF &initDistrib, const VF &gapOpen,
                      const VF &gapExtend, const VVF &emitPairs, const VF &emitSingle, const char *filename);
MultiSequence *DoAlign (MultiSequence *sequence, const ProbabilisticModel &model, VF &initDistrib, VF &gapOpen, VF &gapExtend,
			VVF &emitPairs, VF &emitSingle);
SafeVector<string> ParseParams (int argc, char **argv);
void ReadParameters ();
MultiSequence *ComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,
                                      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                      const ProbabilisticModel &model);
MultiSequence *AlignAlignments (MultiSequence *align1, MultiSequence *align2,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model);
SafeVector<SafeVector<SparseMatrix *> > DoRelaxation (MultiSequence *sequences, 
						      SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, int consgap);
void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior);
void Relax1 (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior);
void Relax_gap (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior);

set<int> GetSubtree (const TreeNode *tree);
void TreeBasedBiPartitioning (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                              const ProbabilisticModel &model, MultiSequence* &alignment,
                              const TreeNode *tree);
void DoIterativeRefinement (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model, MultiSequence* &alignment);
void WriteAnnotation (MultiSequence *alignment,
		      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);
void ComputeAnnotation (MultiSequence *alignment,
          const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);

int ComputeScore (const SafeVector<pair<int, int> > &active, 
		  const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);

/////////////////////////////////////////////////////////////////
// main()
//
// Calls all initialization routines and runs the PROBCONS
// aligner.
/////////////////////////////////////////////////////////////////

int main (int argc, char **argv){

  // print PROBCONS heading
  PrintHeading();
  
  // parse program parameters
  SafeVector<string> sequenceNames = ParseParams (argc, argv);
  ReadParameters();
  PrintParameters ("Using parameter set:", initDistrib, gapOpen, gapExtend, emitPairs, emitSingle, NULL);

  // now, we'll process all the files given as input.  If we are given
  // several filenames as input, then we'll load all of those sequences
  // simultaneously, as long as we're not training.  On the other hand,
  // if we are training, then we'll treat each file as a separate
  // training instance
  
  // if we are training
  if (enableTraining){

    // build new model for aligning
    ProbabilisticModel model (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle, gapSwitchDefault);

    // prepare to average parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] = 0;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] = 0;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] = 0;
    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++)
	     for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] = 0;
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] = 0;
    }
   
    // align each file individually
    for (int i = 0; i < (int) sequenceNames.size(); i++){

      VF thisInitDistrib (NumMatrixTypes);
      VF thisGapOpen (2*NumInsertStates);
      VF thisGapExtend (2*NumInsertStates);
      VVF thisEmitPairs (256, VF (256, 1e-10));
      VF thisEmitSingle (256, 1e-5);
      
      // load sequence file
      MultiSequence *sequences = new MultiSequence(); assert (sequences);
      cerr << "Loading sequence file: " << sequenceNames[i] << endl;
      sequences->LoadMFA (sequenceNames[i], true);

      // align sequences
      DoAlign (sequences, model, thisInitDistrib, thisGapOpen, thisGapExtend, thisEmitPairs, thisEmitSingle);

      // add in contribution of the derived parameters
      for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] += thisInitDistrib[i];
      for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] += thisGapOpen[i];
      for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] += thisGapExtend[i];
      if (enableTrainEmissions){
    	for (int i = 0; i < (int) emitPairs.size(); i++) 
    	  for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] += thisEmitPairs[i][j];
    	for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] += thisEmitSingle[i];
      }

      delete sequences;
    }

    // compute new parameters and print them out
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] /= (int) sequenceNames.size();
    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++) 
      	for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] /= (int) sequenceNames.size();
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] /= sequenceNames.size();
    }
    
    PrintParameters ("Trained parameter set:",
                     initDistrib, gapOpen, gapExtend, emitPairs, emitSingle,
                     parametersOutputFilename.c_str());
  }

  // if we are not training, we must simply want to align some sequences
  else {

    // load all files together
    MultiSequence *sequences = new MultiSequence(); assert (sequences);
    for (int i = 0; i < (int) sequenceNames.size(); i++){
      cerr << "Loading sequence file: " << sequenceNames[i] << endl;
      sequences->LoadMFA (sequenceNames[i], true);
    }

    // do all "pre-training" repetitions first
    for (int ct = 0; ct < numPreTrainingReps; ct++){
      enableTraining = true;

      // build new model for aligning
      ProbabilisticModel model (initDistrib, gapOpen, gapExtend, 
                                emitPairs, emitSingle, gapSwitchDefault);

      // do initial alignments
      DoAlign (sequences, model, initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);

      // print new parameters
      PrintParameters ("Recomputed parameter set:", initDistrib, gapOpen, gapExtend, emitPairs, emitSingle, NULL);

      enableTraining = false;
    }

    // now, we can perform the alignments and write them out
    MultiSequence *alignment = DoAlign (sequences,
                                        ProbabilisticModel (initDistrib, gapOpen, gapExtend,  emitPairs, emitSingle, gapSwitchDefault),
                                        initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);
    
    if (!enableAllPairs){
      if (enableClustalWOutput)
      	alignment->WriteALN (cout);
      else
      	alignment->WriteMFA (cout);
    }
    delete alignment;
    delete sequences;
   
  }
}

/////////////////////////////////////////////////////////////////
// PrintHeading()
//
// Prints heading for PROBCONS program.
/////////////////////////////////////////////////////////////////

void PrintHeading (){
  cerr << endl
       << "PROBCONS version " << VERSION << " - align multiple protein sequences and print to standard output" << endl
       << "Written by Chuong Do" << endl
       << endl;
}

/////////////////////////////////////////////////////////////////
// PrintParameters()
//
// Prints PROBCONS parameters to STDERR.  If a filename is
// specified, then the parameters are also written to the file.
/////////////////////////////////////////////////////////////////

void PrintParameters (const char *message, const VF &initDistrib, const VF &gapOpen, const VF &gapExtend, const VVF &emitPairs, const VF &emitSingle, const char *filename) {
  // print parameters to the screen
  cerr << message << endl
       << "    initDistrib[] = { ";
  for (int i = 0; i < NumMatrixTypes; i++) cerr << setprecision (10) << initDistrib[i] << " ";
  cerr << "}" << endl
       << "        gapOpen[] = { ";
  for (int i = 0; i < NumInsertStates*2; i++) cerr << setprecision (10) << gapOpen[i] << " ";
  cerr << "}" << endl
       << "      gapExtend[] = { ";
  for (int i = 0; i < NumInsertStates*2; i++) cerr << setprecision (10) << gapExtend[i] << " ";
  cerr << "}" << endl
       << endl;

  // if a file name is specified
  if (filename){

    // attempt to open the file for writing
    FILE *file = fopen (filename, "w");
    if (!file){
      cerr << "ERROR: Unable to write parameter file: " << filename << endl;
      exit (1);
    }

    // if successful, then write the parameters to the file
    for (int i = 0; i < NumMatrixTypes; i++) fprintf (file, "%.10f ", initDistrib[i]); fprintf (file, "\n");
    for (int i = 0; i < 2*NumInsertStates; i++) fprintf (file, "%.10f ", gapOpen[i]); fprintf (file, "\n");
    for (int i = 0; i < 2*NumInsertStates; i++) fprintf (file, "%.10f ", gapExtend[i]); fprintf (file, "\n");
    fprintf (file, "%s\n", alphabet.c_str());
    for (int i = 0; i < (int) alphabet.size(); i++){
      for (int j = 0; j <= i; j++)
	fprintf (file, "%.10f ", emitPairs[(unsigned char) alphabet[i]][(unsigned char) alphabet[j]]);
      fprintf (file, "\n");
    }
    for (int i = 0; i < (int) alphabet.size(); i++)
      fprintf (file, "%.10f ", emitSingle[(unsigned char) alphabet[i]]);
    fprintf (file, "\n");
    fclose (file);
  }
}

/////////////////////////////////////////////////////////////////
// DoAlign()
//
// First computes all pairwise posterior probability matrices.
// Then, computes new parameters if training, or a final
// alignment, otherwise.
/////////////////////////////////////////////////////////////////

MultiSequence *DoAlign (MultiSequence *sequences, const ProbabilisticModel &model, VF &initDistrib, VF &gapOpen, VF &gapExtend, VVF &emitPairs, VF &emitSingle){

  assert (sequences);

  const int numSeqs = sequences->GetNumSequences();
  VVF distances (numSeqs, VF (numSeqs, 0));
  SafeVector<SafeVector<SparseMatrix *> > sparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));
  // SafeVector<SafeVector<SparseMatrix *> > untransformedSparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));

  if (enableTraining){
    // prepare to average parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] = 0;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] = 0;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] = 0;
    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++)
        for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] = 0;
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] = 0;
    }
  }

  // skip posterior calculations if we just want to do Viterbi alignments
  if (!enableViterbi){

    // do all pairwise alignments for posterior probability matrices
    for (int a = 0; a < numSeqs-1; a++){
      for (int b = a+1; b < numSeqs; b++){
        Sequence *seq1 = sequences->GetSequence (a);
        Sequence *seq2 = sequences->GetSequence (b);

        // verbose output
        if (enableVerbose)
          cerr << "Computing posterior matrix: (" << a+1 << ") " << seq1->GetHeader() << " vs. "
               << "(" << b+1 << ") " << seq2->GetHeader() << " -- ";

        // compute forward and backward probabilities
        VF *forward = model.ComputeForwardMatrix (seq1, seq2); assert (forward);
        VF *backward = model.ComputeBackwardMatrix (seq1, seq2); assert (backward);

        // if we are training, then we'll simply want to compute the
        // expected counts for each region within the matrix separately;
        // otherwise, we'll need to put all of the regions together and
        // assemble a posterior probability match matrix

        // so, if we're training
        if (enableTraining){
          
          // compute new parameters
          VF thisInitDistrib (NumMatrixTypes);
          VF thisGapOpen (2*NumInsertStates);
          VF thisGapExtend (2*NumInsertStates);
          VVF thisEmitPairs (256, VF (256, 1e-10));
          VF thisEmitSingle (256, 1e-5);
          
          model.ComputeNewParameters (seq1, seq2, *forward, *backward, thisInitDistrib, thisGapOpen, thisGapExtend, thisEmitPairs, thisEmitSingle, enableTrainEmissions);

          // add in contribution of the derived parameters
          for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] += thisInitDistrib[i];
          for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] += thisGapOpen[i];
          for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] += thisGapExtend[i];
          if (enableTrainEmissions){
            for (int i = 0; i < (int) emitPairs.size(); i++) 
              for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] += thisEmitPairs[i][j];
            for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] += thisEmitSingle[i];
          }

          // let us know that we're done.
          if (enableVerbose) cerr << "done." << endl;
        }
        else {

          // compute posterior probability matrix
          VF *posterior = model.ComputePosteriorMatrix (seq1, seq2, *forward, *backward); assert (posterior);

          // compute sparse representations
          sparseMatrices[a][b] = new SparseMatrix (seq1->GetLength(), seq2->GetLength(), *posterior);
          sparseMatrices[b][a] = NULL; 

          // untransformedSparseMatrices[a][b] = new SparseMatrix (seq1->GetLength(), seq2->GetLength(), *posterior);
          // untransformedSparseMatrices[b][a] = NULL;
          
          if (!enableAllPairs){
            // perform the pairwise sequence alignment
            pair<SafeVector<char> *, float> alignment = model.ComputeAlignment (seq1->GetLength(),
        									seq2->GetLength(),
        									*posterior);
            
            // compute "expected accuracy" distance for evolutionary tree computation
            float distance = alignment.second / min (seq1->GetLength(), seq2->GetLength());
            distances[a][b] = distances[b][a] = distance;
            
            if (enableVerbose)
              cerr << setprecision (10) << distance << endl;
            
              delete alignment.first;
          }
          else {
            // let us know that we're done.
            if (enableVerbose) cerr << "done." << endl;
          }
          
          delete posterior;
        }

        delete forward;
        delete backward;
      }
    }
  }

  // now average out parameters derived
  if (enableTraining){

    // compute new parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] /= numSeqs * (numSeqs - 1) / 2;

    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++)
        for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] /= numSeqs * (numSeqs - 1) / 2;
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] /= numSeqs * (numSeqs - 1) / 2;
    }
  }

  // see if we still want to do some alignments
  else {

    if (!enableViterbi){

      // perform the consistency transformation the desired number of times
      for (int r = 0; r < numConsistencyReps; r++) {
        SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices = DoRelaxation (sequences, sparseMatrices, consgap);

        // now replace the old posterior matrices
        for (int i = 0; i < numSeqs; i++) {
          for (int j = 0; j < numSeqs; j++) {
            delete sparseMatrices[i][j];
            sparseMatrices[i][j] = newSparseMatrices[i][j];
          }
        }

      }

    }

    MultiSequence *finalAlignment = NULL;

    if (enableAllPairs){
      for (int a = 0; a < numSeqs-1; a++){
        for (int b = a+1; b < numSeqs; b++){
          Sequence *seq1 = sequences->GetSequence (a);
          Sequence *seq2 = sequences->GetSequence (b);
          
          if (enableVerbose)
            cerr << "Performing pairwise alignment: (" << a+1 << ") " << seq1->GetHeader() << " vs. "
        	 << "(" << b+1 << ") " << seq2->GetHeader() << " -- ";

          
          // perform the pairwise sequence alignment
          pair<SafeVector<char> *, float> alignment;
          if (enableViterbi)
            alignment = model.ComputeViterbiAlignment (seq1, seq2);
          else {

            // build posterior matrix
            VF *posterior = sparseMatrices[a][b]->GetPosterior(); assert (posterior);
            int length = (seq1->GetLength() + 1) * (seq2->GetLength() + 1);
            for (int i = 0; i < length; i++) (*posterior)[i] -= cutoff;

            alignment = model.ComputeAlignment (seq1->GetLength(), seq2->GetLength(), *posterior);
            delete posterior;
          }

          // write pairwise alignments 
          string name = seq1->GetHeader() + "-" + seq2->GetHeader() + (enableClustalWOutput ? ".aln" : ".fasta");
          ofstream outfile (name.c_str());
          
          MultiSequence *result = new MultiSequence();
          result->AddSequence (seq1->AddGaps(alignment.first, 'X'));
          result->AddSequence (seq2->AddGaps(alignment.first, 'Y'));
          if (enableClustalWOutput)
            result->WriteALN (outfile);
          else
            result->WriteMFA (outfile);
          
          outfile.close();
          
          delete alignment.first;
        }
      }
    }
    
    // now if we still need to do a final multiple alignment
    else {
      
      if (enableVerbose)
	       cerr << endl;
      
      // compute the evolutionary tree
      TreeNode *tree = TreeNode::ComputeTree (distances);
      
      // tree->Print (cerr, sequences);
      // cerr << endl;
      
      // make the final alignment
      finalAlignment = ComputeFinalAlignment (tree, sequences, sparseMatrices, model);
      
      // build annotation
      if (enableAnnotation){
        // WriteAnnotation (finalAlignment, sparseMatrices);
        ComputeAnnotation (finalAlignment, sparseMatrices);
        // ComputeAnnotation (finalAlignment, untransformedSparseMatrices);
      }

      delete tree;
    }

    if (!enableViterbi){
      // delete sparse matrices
      for (int a = 0; a < numSeqs-1; a++){
      	for (int b = a+1; b < numSeqs; b++){
      	  delete sparseMatrices[a][b];
      	  delete sparseMatrices[b][a];
      	}
      }
    }

    return finalAlignment;
  }

  return NULL;
}

/////////////////////////////////////////////////////////////////
// GetInteger()
//
// Attempts to parse an integer from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool GetInteger (char *data, int *val){
  char *endPtr;
  long int retVal;

  assert (val);

  errno = 0;
  retVal = strtol (data, &endPtr, 0);
  if (retVal == 0 && (errno != 0 || data == endPtr)) return false;
  if (errno != 0 && (retVal == LONG_MAX || retVal == LONG_MIN)) return false;
  if (retVal < (long) INT_MIN || retVal > (long) INT_MAX) return false;
  *val = (int) retVal;
  return true;
}

/////////////////////////////////////////////////////////////////
// GetFloat()
//
// Attempts to parse a float from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool GetFloat (char *data, float *val){
  char *endPtr;
  double retVal;

  assert (val);

  errno = 0;
  retVal = strtod (data, &endPtr);
  if (retVal == 0 && (errno != 0 || data == endPtr)) return false;
  if (errno != 0 && (retVal >= 1000000.0 || retVal <= -1000000.0)) return false;
  *val = (float) retVal;
  return true;
}

/////////////////////////////////////////////////////////////////
// ParseParams()
//
// Parse all command-line options.
/////////////////////////////////////////////////////////////////

SafeVector<string> ParseParams (int argc, char **argv){

  if (argc < 2){

    cerr << "PROBCONS comes with ABSOLUTELY NO WARRANTY.  This is free software, and" << endl
         << "you are welcome to redistribute it under certain conditions.  See the" << endl
         << "file COPYING.txt for details." << endl
         << endl
         << "Usage:" << endl
         << "       probcons [OPTION]... [MFAFILE]..." << endl
         << endl
         << "Description:" << endl
         << "       Align sequences in MFAFILE(s) and print result to standard output" << endl
         << endl
         << "       -clustalw" << endl
         << "              use CLUSTALW output format instead of MFA" << endl
         << endl
         << "       -c, --consistency REPS" << endl
         << "              use " << MIN_CONSISTENCY_REPS << " <= REPS <= " << MAX_CONSISTENCY_REPS
         << " (default: " << numConsistencyReps << ") passes of consistency transformation" << endl
         << endl
         << "       -ir, --iterative-refinement REPS" << endl
         << "              use " << MIN_ITERATIVE_REFINEMENT_REPS << " <= REPS <= " << MAX_ITERATIVE_REFINEMENT_REPS
         << " (default: " << numIterativeRefinementReps << ") passes of iterative-refinement" << endl
         << endl
	 << "       -pre, --pre-training REPS" << endl
	 << "              use " << MIN_PRETRAINING_REPS << " <= REPS <= " << MAX_PRETRAINING_REPS
	 << " (default: " << numPreTrainingReps << ") rounds of pretraining" << endl
	 << endl
	 << "       -pairs" << endl
         << "              generate all-pairs pairwise alignments" << endl
         << endl
	 << "       -viterbi" << endl
	 << "              use Viterbi algorithm to generate all pairs (automatically enables -pairs)" << endl
	 << endl
         << "       -v, --verbose" << endl
         << "              report progress while aligning (default: " << (enableVerbose ? "on" : "off") << ")" << endl
         << endl
         << "       -annot FILENAME" << endl
         << "              write annotation for multiple alignment to FILENAME" << endl
         << endl
         << "       -t, --train FILENAME" << endl
         << "              compute EM transition probabilities, store in FILENAME (default: "
         << parametersOutputFilename << ")" << endl
         << endl
         << "       -e, --emissions" << endl
         << "              also reestimate emission probabilities (default: "
         << (enableTrainEmissions ? "on" : "off") << ")" << endl
         << endl
	 << "       -p, --paramfile FILENAME" << endl
	 << "              read parameters from FILENAME (default: "
	 << parametersInputFilename << ")" << endl
	 << endl
	 << "       -a, --alignment-order" << endl
	 << "              print sequences in alignment order rather than input order (default: "
	 << (enableAlignOrder ? "on" : "off") << ")" << endl
	 << endl;
    //     	 << "       -go, --gap-open VALUE" << endl
    //     	 << "              gap opening penalty of VALUE <= 0 (default: " << gapOpenPenalty << ")" << endl
    //	 << endl
    //	 << "       -ge, --gap-extension VALUE" << endl
    //	 << "              gap extension penalty of VALUE <= 0 (default: " << gapContinuePenalty << ")" << endl
    //	 << endl
    //         << "       -co, --cutoff CUTOFF" << endl
    //         << "              subtract 0 <= CUTOFF <= 1 (default: " << cutoff << ") from all posterior values before final alignment" << endl
    //         << endl
    
    exit (1);
  }

  SafeVector<string> sequenceNames;
  int tempInt;
  float tempFloat;

  for (int i = 1; i < argc; i++){
    if (argv[i][0] == '-'){

      // training
      if (!strcmp (argv[i], "-t") || !strcmp (argv[i], "--train")){
        enableTraining = true;
        if (i < argc - 1)
          parametersOutputFilename = string (argv[++i]);
        else {
          cerr << "ERROR: Filename expected for option " << argv[i] << endl;
          exit (1);
        }
      }
      
      // emission training
      else if (!strcmp (argv[i], "-e") || !strcmp (argv[i], "--emissions")){
        enableTrainEmissions = true;
      }

      // parameter file
      else if (!strcmp (argv[i], "-p") || !strcmp (argv[i], "--paramfile")){
        if (i < argc - 1)
          parametersInputFilename = string (argv[++i]);
        else {
          cerr << "ERROR: Filename expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // number of consistency transformations
      else if (!strcmp (argv[i], "-c") || !strcmp (argv[i], "--consistency")){
        if (i < argc - 1){
          if (!GetInteger (argv[++i], &tempInt)){
            cerr << "ERROR: Invalid integer following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempInt < MIN_CONSISTENCY_REPS || tempInt > MAX_CONSISTENCY_REPS){
              cerr << "ERROR: For option " << argv[i-1] << ", integer must be between "
                   << MIN_CONSISTENCY_REPS << " and " << MAX_CONSISTENCY_REPS << "." << endl;
              exit (1);
            }
            else
              numConsistencyReps = tempInt;
          }
        }
        else {
          cerr << "ERROR: Integer expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // number of randomized partitioning iterative refinement passes
      else if (!strcmp (argv[i], "-ir") || !strcmp (argv[i], "--iterative-refinement")){
        if (i < argc - 1){
          if (!GetInteger (argv[++i], &tempInt)){
            cerr << "ERROR: Invalid integer following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempInt < MIN_ITERATIVE_REFINEMENT_REPS || tempInt > MAX_ITERATIVE_REFINEMENT_REPS){
              cerr << "ERROR: For option " << argv[i-1] << ", integer must be between "
                   << MIN_ITERATIVE_REFINEMENT_REPS << " and " << MAX_ITERATIVE_REFINEMENT_REPS << "." << endl;
              exit (1);
            }
            else
              numIterativeRefinementReps = tempInt;
          }
        }
        else {
          cerr << "ERROR: Integer expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // number of EM pre-training rounds
      else if (!strcmp (argv[i], "-pre") || !strcmp (argv[i], "--pre-training")){
        if (i < argc - 1){
          if (!GetInteger (argv[++i], &tempInt)){
            cerr << "ERROR: Invalid integer following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempInt < MIN_PRETRAINING_REPS || tempInt > MAX_PRETRAINING_REPS){
              cerr << "ERROR: For option " << argv[i-1] << ", integer must be between "
                   << MIN_PRETRAINING_REPS << " and " << MAX_PRETRAINING_REPS << "." << endl;
              exit (1);
            }
            else
              numPreTrainingReps = tempInt;
          }
        }
        else {
          cerr << "ERROR: Integer expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // gap open penalty
      else if (!strcmp (argv[i], "-go") || !strcmp (argv[i], "--gap-open")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat > 0){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must not be positive." << endl;
              exit (1);
            }
            else
              gapOpenPenalty = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // gap extension penalty
      else if (!strcmp (argv[i], "-ge") || !strcmp (argv[i], "--gap-extension")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat > 0){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must not be positive." << endl;
              exit (1);
            }
            else
              gapContinuePenalty = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // all-pairs pairwise alignments
      else if (!strcmp (argv[i], "-pairs")){
        enableAllPairs = true;
      }

      // all-pairs pairwise Viterbi alignments
      else if (!strcmp (argv[i], "-viterbi")){
        enableAllPairs = true;
	enableViterbi = true;
      }

      // annotation files
      else if (!strcmp (argv[i], "-annot")){
        enableAnnotation = true;
        if (i < argc - 1)
	  annotationFilename = argv[++i];
        else {
          cerr << "ERROR: FILENAME expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // clustalw output format
      else if (!strcmp (argv[i], "-clustalw")){
	enableClustalWOutput = true;
      }

      // cutoff
      else if (!strcmp (argv[i], "-co") || !strcmp (argv[i], "--cutoff")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat < 0 || tempFloat > 1){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must be between 0 and 1." << endl;
              exit (1);
            }
            else
              cutoff = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // verbose reporting
      else if (!strcmp (argv[i], "-v") || !strcmp (argv[i], "--verbose")){
        enableVerbose = true;
      }

      // alignment order
      else if (!strcmp (argv[i], "-a") || !strcmp (argv[i], "--alignment-order")){
	enableAlignOrder = true;
      }

      // bad arguments
      else {
        cerr << "ERROR: Unrecognized option: " << argv[i] << endl;
        exit (1);
      }
    }
    else {
      sequenceNames.push_back (string (argv[i]));
    }
  }

  if (enableTrainEmissions && !enableTraining){
    cerr << "ERROR: Training emissions (-e) requires training (-t)" << endl;
    exit (1);
  }

  return sequenceNames;
}

/////////////////////////////////////////////////////////////////
// ReadParameters()
//
// Read initial distribution, transition, and emission
// parameters from a file.
/////////////////////////////////////////////////////////////////

void ReadParameters (){

  ifstream data;

  emitPairs = VVF (256, VF (256, 1e-10));
  emitSingle = VF (256, 1e-5);

  // read initial state distribution and transition parameters
  if (parametersInputFilename == string ("")){
    if (NumInsertStates == 1){
      for (int i = 0; i < NumMatrixTypes; i++) initDistrib[i] = initDistrib1Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapOpen[i] = gapOpen1Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapExtend[i] = gapExtend1Default[i];
    }
    else if (NumInsertStates == 2){
      for (int i = 0; i < NumMatrixTypes; i++) initDistrib[i] = initDistrib2Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapOpen[i] = gapOpen2Default[i];
      for (int i = 0; i < 2*NumInsertStates; i++) gapExtend[i] = gapExtend2Default[i];
    }
    else {
      cerr << "ERROR: No default initial distribution/parameter settings exist" << endl
           << "       for " << NumInsertStates << " pairs of insert states.  Use --paramfile." << endl;
      exit (1);
    }

    alphabet = alphabetDefault;

    for (int i = 0; i < (int) alphabet.length(); i++){
      emitSingle[(unsigned char) tolower(alphabet[i])] = emitSingleDefault[i];
      emitSingle[(unsigned char) toupper(alphabet[i])] = emitSingleDefault[i];
      for (int j = 0; j <= i; j++){
        emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(alphabet[j])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(alphabet[i])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(alphabet[i])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(alphabet[i])] = emitPairsDefault[i][j];
        emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(alphabet[i])] = emitPairsDefault[i][j];
      }
    }
  }
  else {
    data.open (parametersInputFilename.c_str());
    if (data.fail()){
      cerr << "ERROR: Unable to read parameter file: " << parametersInputFilename << endl;
      exit (1);
    }
    
    string line[3];
    for (int i = 0; i < 3; i++){
      if (!getline (data, line[i])){
        cerr << "ERROR: Unable to read transition parameters from parameter file: " << parametersInputFilename << endl;
        exit (1);
      }
    }
    istringstream data2;
    data2.clear(); data2.str (line[0]); for (int i = 0; i < NumMatrixTypes; i++) data2 >> initDistrib[i];
    data2.clear(); data2.str (line[1]); for (int i = 0; i < 2*NumInsertStates; i++) data2 >> gapOpen[i];
    data2.clear(); data2.str (line[2]); for (int i = 0; i < 2*NumInsertStates; i++) data2 >> gapExtend[i];

    if (!getline (data, line[0])){
      cerr << "ERROR: Unable to read alphabet from scoring matrix file: " << parametersInputFilename << endl;
      exit (1);
    }
    
    // read alphabet as concatenation of all characters on alphabet line
    alphabet = "";
    string token;
    data2.clear(); data2.str (line[0]); while (data2 >> token) alphabet += token;

    for (int i = 0; i < (int) alphabet.size(); i++){
      for (int j = 0; j <= i; j++){
        float val;
        data >> val;
        emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(alphabet[j])] = val;
        emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(alphabet[j])] = val;
        emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(alphabet[j])] = val;
        emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(alphabet[j])] = val;
        emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(alphabet[i])] = val;
        emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(alphabet[i])] = val;
        emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(alphabet[i])] = val;
        emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(alphabet[i])] = val;
      }
    }

    for (int i = 0; i < (int) alphabet.size(); i++){
      float val;
      data >> val;
      emitSingle[(unsigned char) tolower(alphabet[i])] = val;
      emitSingle[(unsigned char) toupper(alphabet[i])] = val;
    }
    data.close();
  }
}

/////////////////////////////////////////////////////////////////
// ProcessTree()
//
// Process the tree recursively.  Returns the aligned sequences
// corresponding to a node or leaf of the tree.
/////////////////////////////////////////////////////////////////

MultiSequence *ProcessTree (const TreeNode *tree, MultiSequence *sequences,const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,const ProbabilisticModel &model){
  MultiSequence *result;

  // check if this is a node of the alignment tree
  if (tree->GetSequenceLabel() == -1){
    MultiSequence *alignLeft = ProcessTree (tree->GetLeftChild(), sequences, sparseMatrices, model);
    MultiSequence *alignRight = ProcessTree (tree->GetRightChild(), sequences, sparseMatrices, model);

    assert (alignLeft);
    assert (alignRight);

    result = AlignAlignments (alignLeft, alignRight, sparseMatrices, model);
    assert (result);

    delete alignLeft;
    delete alignRight;
  }

  // otherwise, this is a leaf of the alignment tree
  else {
    result = new MultiSequence(); assert (result);
    result->AddSequence (sequences->GetSequence(tree->GetSequenceLabel())->Clone());
  }

  return result;
}

/////////////////////////////////////////////////////////////////
// ComputeFinalAlignment()
//
// Compute the final alignment by calling ProcessTree(), then
// performing iterative refinement as needed.
/////////////////////////////////////////////////////////////////

MultiSequence *ComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,const ProbabilisticModel &model){

  MultiSequence *alignment = ProcessTree (tree, sequences, sparseMatrices, model);

  SafeVector<int> oldOrdering;
  if (enableAlignOrder){
    for (int i = 0; i < alignment->GetNumSequences(); i++)
      oldOrdering.push_back (alignment->GetSequence(i)->GetSortLabel());
    alignment->SaveOrdering();
    enableAlignOrder = false;
  }

  // tree-based refinement
  // TreeBasedBiPartitioning (sparseMatrices, model, alignment, tree);

  // iterative refinement
  for (int i = 0; i < numIterativeRefinementReps; i++)
    DoIterativeRefinement (sparseMatrices, model, alignment);

  // cerr << endl;

  if (oldOrdering.size() > 0){
    for (int i = 0; i < (int) oldOrdering.size(); i++){
      alignment->GetSequence(i)->SetSortLabel(oldOrdering[i]);
    }
  }

  // return final alignment
  return alignment;
}

/////////////////////////////////////////////////////////////////
// AlignAlignments()
//
// Returns the alignment of two MultiSequence objects.
/////////////////////////////////////////////////////////////////

MultiSequence *AlignAlignments (MultiSequence *align1, MultiSequence *align2,const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,const ProbabilisticModel &model){

  // print some info about the alignment
  if (enableVerbose){
    for (int i = 0; i < align1->GetNumSequences(); i++)
      cerr << ((i==0) ? "[" : ",") << align1->GetSequence(i)->GetLabel();
    cerr << "] vs. ";
    for (int i = 0; i < align2->GetNumSequences(); i++)
      cerr << ((i==0) ? "[" : ",") << align2->GetSequence(i)->GetLabel();
    cerr << "]: ";
  }

  VF *posterior = model.BuildPosterior (align1, align2, sparseMatrices, cutoff);
  pair<SafeVector<char> *, float> alignment;

  // choose the alignment routine depending on the "cosmetic" gap penalties used
  if (gapOpenPenalty == 0 && gapContinuePenalty == 0)
    alignment = model.ComputeAlignment (align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), *posterior);
  else
    alignment = model.ComputeAlignmentWithGapPenalties (align1, align2,
                                                        *posterior, align1->GetNumSequences(), align2->GetNumSequences(),
                                                        gapOpenPenalty, gapContinuePenalty);

  delete posterior;

  if (enableVerbose){

    // compute total length of sequences
    int totLength = 0;
    for (int i = 0; i < align1->GetNumSequences(); i++)
      for (int j = 0; j < align2->GetNumSequences(); j++)
        totLength += min (align1->GetSequence(i)->GetLength(), align2->GetSequence(j)->GetLength());

    // give an "accuracy" measure for the alignment
    cerr << alignment.second / totLength << endl;
  }

  // now build final alignment
  MultiSequence *result = new MultiSequence();
  for (int i = 0; i < align1->GetNumSequences(); i++)
    result->AddSequence (align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
  for (int i = 0; i < align2->GetNumSequences(); i++)
    result->AddSequence (align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));
  if (!enableAlignOrder)
    result->SortByLabel();

  // free temporary alignment
  delete alignment.first;

  return result;
}

/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the consistency transformation.  The
// formula used is:
//                     1
//    P'(x[i]-y[j]) = ---  sum   sum P(x[i]-z[k]) P(z[k]-y[j])
//                    |S| z in S  k
//
// where S = {x, y, all other sequences...}
//
/////////////////////////////////////////////////////////////////

SafeVector<SafeVector<SparseMatrix *> > DoRelaxation (MultiSequence *sequences, SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, int consgap){
  const int numSeqs = sequences->GetNumSequences();

  SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));

  // for every pair of sequences
  for (int i = 0; i < numSeqs; i++){
    for (int j = i+1; j < numSeqs; j++){

      // cerr << "Relaxing all pairs --> i:" << i << " j:" << j << endl;

      Sequence *seq1 = sequences->GetSequence (i);
      Sequence *seq2 = sequences->GetSequence (j);

      if (enableVerbose)
        cerr << "Relaxing (" << i+1 << ") " << seq1->GetHeader() << " vs. "
             << "(" << j+1 << ") " << seq2->GetHeader() << ": ";

      // get the original posterior matrix
      VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior(); assert (posteriorPtr);
      VF &posterior = *posteriorPtr;

      const int seq1Length = seq1->GetLength();
      const int seq2Length = seq2->GetLength();

      // contribution from the summation where z = x and z = y
      for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] += posterior[k];

      if (enableVerbose)
        cerr << sparseMatrices[i][j]->GetNumCells() << " --> ";

      // contribution from all other sequences
      for (int k = 0; k < numSeqs; k++) if (k != i && k != j){

        if (k < i) {
          // cerr << "Relax1 k:" << k << " i:" << i << " j:" << j << endl;
          SparseMatrix *temp = sparseMatrices[k][i]->ComputeTranspose();
          if (consgap) {
            Relax_gap (temp, sparseMatrices[k][j], posterior);
          }
          else {
            Relax (temp, sparseMatrices[k][j], posterior);
          }
          // Relax1 (sparseMatrices[k][i], sparseMatrices[k][j], posterior);
          delete temp;
        } else if (k > i && k < j){
          // cerr << "Relax k:" << k << " i:" << i << " j:" << j << endl;
          if (consgap) {
            Relax_gap (sparseMatrices[i][k], sparseMatrices[k][j], posterior);
          }
          else {
            Relax (sparseMatrices[i][k], sparseMatrices[k][j], posterior);
          }
        } else {
          // cerr << "Relax temp k:" << k << " i:" << i << " j:" << j << endl;
          SparseMatrix *temp = sparseMatrices[j][k]->ComputeTranspose();
          if (consgap) {
            Relax_gap (sparseMatrices[i][k], temp, posterior);
          } else {
            Relax (sparseMatrices[i][k], temp, posterior);
          }
          delete temp;
        }
      }

      // now renormalization
      for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] /= numSeqs;

      // mask out positions not originally in the posterior matrix
      SparseMatrix *matXY = sparseMatrices[i][j];
      
      for (int y = 0; y <= seq2Length; y++) posterior[y] = 0;
      
      for (int x = 1; x <= seq1Length; x++){
        SafeVector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
        SafeVector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
        VF::iterator base = posterior.begin() + x * (seq2Length + 1);
        int curr = 0;

        while (XYptr != XYend){
          // zero out all cells until the first filled column
          while (curr < XYptr->first){
            base[curr] = 0;
            curr++;
          }

          // now, skip over this column
          curr++;
          ++XYptr;
        }

        // zero out cells after last column
        while (curr <= seq2Length){
          base[curr] = 0;
          curr++;
        }
      }

      // save the new posterior matrix
      newSparseMatrices[i][j] = new SparseMatrix (seq1->GetLength(), seq2->GetLength(), posterior);
      newSparseMatrices[j][i] = NULL;

      if (enableVerbose)
        cerr << newSparseMatrices[i][j]->GetNumCells() << " -- ";

      delete posteriorPtr;

      if (enableVerbose)
        cerr << "done." << endl;
    }
  }
  
  return newSparseMatrices;
}


void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior){

  assert (matXZ);
  assert (matZY);

  int lengthX = matXZ->GetSeq1Length();
  int lengthY = matZY->GetSeq2Length();
  assert (matXZ->GetSeq2Length() == matZY->GetSeq1Length());

  // for every x[i]
  for (int i = 1; i <= lengthX; i++){
    SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
    SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

    VF::iterator base = posterior.begin() + i * (lengthY + 1);

    // iterate through all x[i]-z[k]
    while (XZptr != XZend){
      SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
      SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
      const float XZval = XZptr->second;

      // iterate through all z[k]-y[j]
      while (ZYptr != ZYend){
        base[ZYptr->first] += XZval * ZYptr->second;
        ZYptr++;
      }
      XZptr++;
    }
  }
}


/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax_gap (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior){

  assert (matXZ);
  assert (matZY);

  int lengthX = matXZ->GetSeq1Length();
  int lengthY = matZY->GetSeq2Length();
  assert (matXZ->GetSeq2Length() == matZY->GetSeq1Length());

  float * p_gapX = new float [lengthX+1];
  float * p_gapY = new float [lengthY+1];

  // for every x[i]
  for (int i = 1; i <= lengthX; i++){
    SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
    SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

    VF::iterator base = posterior.begin() + i * (lengthY + 1);

    float totprobXiZ=0; //total probability that X[i] is aligned to any base in Z

    // iterate through all x[i]-z[k]
    while (XZptr != XZend){
      SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
      SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
      const float XZval = XZptr->second;

      totprobXiZ+=XZval;

      // iterate through all z[k]-y[j]
      while (ZYptr != ZYend){
        base[ZYptr->first] += XZval * ZYptr->second;
        ZYptr++;
      }
      XZptr++;
    }
    p_gapX[i]=1-totprobXiZ; //probability that Xi is aligned to a gap in Z
  }

  SparseMatrix *matYZ = matZY->ComputeTranspose();
  for (int j = 1; j <= lengthY; j++){
    float totprobYjZ=0;
    SafeVector<PIF>::iterator YZptr = matYZ->GetRowPtr(j);
    SafeVector<PIF>::iterator YZend = YZptr + matYZ->GetRowSize(j);
    while (YZptr != YZend){
      totprobYjZ+=YZptr->second;
      YZptr++;
    }
    p_gapY[j]=1-totprobYjZ;
  }

  //correct for the probability that both X[i] and Y[j] are aligned to a gap in Z
  for (int i = 1; i <= lengthX; i++){
    VF::iterator base = posterior.begin() + i * (lengthY + 1);
    for (int j = 1; j <= lengthY; j++){
      base[j]+= p_gapX[i] * p_gapY[j];
    }
  }

  delete p_gapX;
  delete p_gapY;
  delete matYZ;
}


//TODO: check this Relax1 function and apply same correction here!


/////////////////////////////////////////////////////////////////
// Relax1()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax1 (SparseMatrix *matZX, SparseMatrix *matZY, VF &posterior){

  assert (matZX);
  assert (matZY);

  int lengthZ = matZX->GetSeq1Length();
  int lengthY = matZY->GetSeq2Length();

  // for every z[k]
  for (int k = 1; k <= lengthZ; k++){
    SafeVector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
    SafeVector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

    // iterate through all z[k]-x[i]
    while (ZXptr != ZXend){
      SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
      SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
      const float ZXval = ZXptr->second;
      VF::iterator base = posterior.begin() + ZXptr->first * (lengthY + 1);

      // iterate through all z[k]-y[j]
      while (ZYptr != ZYend){
        base[ZYptr->first] += ZXval * ZYptr->second;
        ZYptr++;
      }
      ZXptr++;
    }
  }
  
}

/////////////////////////////////////////////////////////////////
// GetSubtree
//
// Returns set containing all leaf labels of the current subtree.
/////////////////////////////////////////////////////////////////

set<int> GetSubtree (const TreeNode *tree){
  set<int> s;

  if (tree->GetSequenceLabel() == -1){
    s = GetSubtree (tree->GetLeftChild());
    set<int> t = GetSubtree (tree->GetRightChild());

    for (set<int>::iterator iter = t.begin(); iter != t.end(); ++iter)
      s.insert (*iter);
  }
  else {
    s.insert (tree->GetSequenceLabel());
  }

  return s;
}

/////////////////////////////////////////////////////////////////
// TreeBasedBiPartitioning
//
// Uses the iterative refinement scheme from MUSCLE.
/////////////////////////////////////////////////////////////////

void TreeBasedBiPartitioning (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,const ProbabilisticModel &model, MultiSequence* &alignment,const TreeNode *tree){
  // check if this is a node of the alignment tree
  if (tree->GetSequenceLabel() == -1){
    TreeBasedBiPartitioning (sparseMatrices, model, alignment, tree->GetLeftChild());
    TreeBasedBiPartitioning (sparseMatrices, model, alignment, tree->GetRightChild());

    set<int> leftSubtree = GetSubtree (tree->GetLeftChild());
    set<int> rightSubtree = GetSubtree (tree->GetRightChild());
    set<int> leftSubtreeComplement, rightSubtreeComplement;

    // calculate complement of each subtree
    for (int i = 0; i < alignment->GetNumSequences(); i++){
      if (leftSubtree.find(i) == leftSubtree.end()) leftSubtreeComplement.insert (i);
      if (rightSubtree.find(i) == rightSubtree.end()) rightSubtreeComplement.insert (i);
    }

    // perform realignments for edge to left child
    if (!leftSubtree.empty() && !leftSubtreeComplement.empty()){
      MultiSequence *groupOneSeqs = alignment->Project (leftSubtree); assert (groupOneSeqs);
      MultiSequence *groupTwoSeqs = alignment->Project (leftSubtreeComplement); assert (groupTwoSeqs);
      delete alignment;
      alignment = AlignAlignments (groupOneSeqs, groupTwoSeqs, sparseMatrices, model);
    }

    // perform realignments for edge to right child
    if (!rightSubtree.empty() && !rightSubtreeComplement.empty()){
      MultiSequence *groupOneSeqs = alignment->Project (rightSubtree); assert (groupOneSeqs);
      MultiSequence *groupTwoSeqs = alignment->Project (rightSubtreeComplement); assert (groupTwoSeqs);
      delete alignment;
      alignment = AlignAlignments (groupOneSeqs, groupTwoSeqs, sparseMatrices, model);
    }
  }
}

/////////////////////////////////////////////////////////////////
// DoIterativeRefinement()
//
// Performs a single round of randomized partionining iterative
// refinement.
/////////////////////////////////////////////////////////////////

void DoIterativeRefinement (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,const ProbabilisticModel &model, MultiSequence* &alignment){
  set<int> groupOne, groupTwo;

  // create two separate groups
  for (int i = 0; i < alignment->GetNumSequences(); i++){
    if (rand() % 2)
      groupOne.insert (i);
    else
      groupTwo.insert (i);
  }

  if (groupOne.empty() || groupTwo.empty()) return;

  // project into the two groups
  MultiSequence *groupOneSeqs = alignment->Project (groupOne); assert (groupOneSeqs);
  MultiSequence *groupTwoSeqs = alignment->Project (groupTwo); assert (groupTwoSeqs);
  delete alignment;

  // realign
  alignment = AlignAlignments (groupOneSeqs, groupTwoSeqs, sparseMatrices, model);

  delete groupOneSeqs;
  delete groupTwoSeqs;
}

/////////////////////////////////////////////////////////////////
// WriteAnnotation()
//
// Computes annotation for multiple alignment and write values
// to a file.
/////////////////////////////////////////////////////////////////

void WriteAnnotation (MultiSequence *alignment, const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices){
  ofstream outfile (annotationFilename.c_str());
  
  if (outfile.fail()){
    cerr << "ERROR: Unable to write annotation file." << endl;
    exit (1);
  }

  const int alignLength = alignment->GetSequence(0)->GetLength();
  const int numSeqs = alignment->GetNumSequences();
  
  SafeVector<int> position (numSeqs, 0);
  SafeVector<SafeVector<char>::iterator> seqs (numSeqs);
  for (int i = 0; i < numSeqs; i++) seqs[i] = alignment->GetSequence(i)->GetDataPtr();
  SafeVector<pair<int,int> > active;
  active.reserve (numSeqs);

  SafeVector<int> lab;
  for (int i = 0; i < numSeqs; i++) lab.push_back(alignment->GetSequence(i)->GetSortLabel());
  
  // for every column
  for (int i = 1; i <= alignLength; i++){
    
    // find all aligned residues in this particular column
    active.clear();
    for (int j = 0; j < numSeqs; j++){
      if (seqs[j][i] != '-'){
	active.push_back (make_pair(lab[j], ++position[j]));
      }
    }
    
    sort (active.begin(), active.end());
    outfile << setw(4) << ComputeScore (active, sparseMatrices) << endl;
  }
  
  outfile.close();
}


/////////////////////////////////////////////////////////////////
// ComputeAnnotation()
//
// Computes annotation for multiple alignment and returns values
// in a vector.
/////////////////////////////////////////////////////////////////

void ComputeAnnotation (MultiSequence *alignment, const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices){

  const int alignLength = alignment->GetSequence(0)->GetLength();
  const int numSeqs = alignment->GetNumSequences();
  
  SafeVector<int> annot;
  SafeVector<int> position (numSeqs, 0);
  SafeVector<SafeVector<char>::iterator> seqs (numSeqs);
  for (int i = 0; i < numSeqs; i++) seqs[i] = alignment->GetSequence(i)->GetDataPtr();
  SafeVector<pair<int,int> > active;
  active.reserve (numSeqs);

  SafeVector<int> lab;
  for (int i = 0; i < numSeqs; i++) lab.push_back(alignment->GetSequence(i)->GetSortLabel());
  
  // for every column
  for (int i = 1; i <= alignLength; i++){
    
    // find all aligned residues in this particular column
    active.clear();
    for (int j = 0; j < numSeqs; j++){
      if (seqs[j][i] != '-'){
        active.push_back (make_pair(lab[j], ++position[j]));
      }
    }
    sort (active.begin(), active.end());
    alignment->conf.push_back(ComputeScore(active, sparseMatrices));
  }
}

/////////////////////////////////////////////////////////////////
// ComputeScore()
//
// Computes the annotation score for a particular column.
/////////////////////////////////////////////////////////////////

int ComputeScore (const SafeVector<pair<int, int> > &active, const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices){

  if (active.size() <= 1) return 0;
  
  // ALTERNATIVE #1: Compute the average alignment score.
  // cerr << "active size" << active.size() << endl;
  float val = 0;
  for (int i = 0; i < (int) active.size(); i++){
    for (int j = i+1; j < (int) active.size(); j++){
      // cerr << "i " << i << " - j " << j << endl;
      // cerr << "i first " << active[i].first << " - j first " << active[j].first << endl;
      // cerr << "i second " << active[i].second << " - j second " << active[j].second << endl;
      val += sparseMatrices[active[i].first][active[j].first]->GetValue(active[i].second, active[j].second);
      // cerr << "value at sparseMatrices[active[i].first][active[j].first] - " << sparseMatrices[active[i].first][active[j].first]->GetValue(active[i].second, active[j].second) << endl;
    }
  }

  // cerr << "total " << val << endl;
  // cerr << "score " << (200 * val / ((int) active.size() * ((int) active.size() - 1))) << endl;

  return (int) (200 * val / ((int) active.size() * ((int) active.size() - 1)));
}



//PYTHON LIB!

static PyObject *bpalign(Probcons *self, PyObject *args, PyObject *keywds)
{
    PyObject * seqlist;
    
    const char *kwlist[] = {"name_seq_list","consistency","refinement","pretraining","consgap",NULL};

    numConsistencyReps = 0;
    numIterativeRefinementReps = 0;
    numPreTrainingReps = 0;
    consgap = 0;

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iiii", const_cast<char **>(kwlist), &seqlist, &numConsistencyReps, &numIterativeRefinementReps, &numPreTrainingReps, &consgap))
        return NULL;
    
    // if (!PyArg_ParseTuple(args, "O", &seqlist))
        // return NULL;

    PyObject *iter;
    PyObject *tup;

    iter=PyObject_GetIter(seqlist);

    char * name;
    char * seq;
    int length,i,seqlabel=0;
    
    // parametersInputFilename="nw.txt";
    
    ReadParameters(); //Takes default hardcoded parameters

    // PrintParameters ("Using parameter set:", initDistrib, gapOpen, gapExtend, emitPairs, emitSingle, NULL);

    MultiSequence *sequences = new MultiSequence(); assert (sequences);

    // sequences->LoadMFA("small.fasta");

    while ((tup=PyIter_Next(iter))){
        PyArg_ParseTuple(tup,"ss#",&name, &seq, &length);
        SafeVector<char> *data;
        data = new SafeVector<char>;
        data->push_back ('@');

        for (i=0; i<length; i++){
          data->push_back (seq[i]);
        }

        Sequence *sequence = new Sequence(data, name, length, seqlabel, seqlabel); //int sequenceLabel, int inputLabel);
        seqlabel++;
        sequences->AddSequence (sequence);
    }

    enableAnnotation=true;

    //do the consistency based multiple sequence alignment
    MultiSequence *alignment = DoAlign (sequences,
                                        ProbabilisticModel (initDistrib, gapOpen, gapExtend,  emitPairs, emitSingle, gapSwitchDefault),
                                        initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);
    
    // alignment->WriteALN (cout);

    Sequence * aln;
    int j=0;
    int cols=0;
    const int numSeqs = alignment->GetNumSequences();

    PyObject *res;
    res=PyList_New(alignment->GetNumSequences());

    //obtain aligned sequence
    for (i=0; i<numSeqs; i++){
      aln=alignment->GetSequence(i);

      SafeVector<char>::iterator data;
      string header;
      string alnseq;

      data = aln->GetDataPtr();
      header = aln->GetHeader();
      cols=aln->GetLength();

      for (j=1 ; j<=cols; j++){
        alnseq.push_back(data[j]);
      }

      PyObject *tup=Py_BuildValue("(s,s)",header.c_str(),alnseq.c_str());

      PyList_SetItem(res,i,tup);
    }
    
    PyObject *conf=PyList_New(cols);

    SafeVector<int>::iterator confiter;
    confiter=alignment->conf.begin();

    for (i=0; i<cols; i++){
      PyObject *confi = Py_BuildValue("i",confiter[i]);
      PyList_SetItem(conf,i,confi);
    }

    PyObject *ret=Py_BuildValue("(O,O)",res,conf);

    delete alignment;
    delete sequences;

    return ret;
}


static PyMethodDef probcons_methods[] = {
    { "align", (PyCFunction) bpalign, METH_VARARGS|METH_KEYWORDS },
    { NULL, NULL }
};

static int
probcons_init(Probcons *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

static PyObject *
probcons_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Probcons *self;

    self = (Probcons *)type->tp_alloc(type, 0);

    if (self!=NULL) {
        //pre-init here...
    }
    
    return (PyObject *)self;
}

static PyGetSetDef probcons_getseters[] = {
    {NULL}
};

static void
probcons_dealloc(Probcons *self)
{
    // fprintf(stderr,"Dealloc Probcons object.\n");
}

static PyTypeObject ProbconsType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "probcons",            /*tp_name*/
    sizeof(Probcons),       /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)probcons_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Probcons",            /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    probcons_methods,            /* tp_methods */
    0,                         /* tp_members */
    probcons_getseters,          /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)probcons_init,     /* tp_init */
    0,                         /* tp_alloc */
    probcons_new,                         /* tp_new */
};

Probcons* newProbcons()
{
    return (Probcons *) PyObject_CallObject((PyObject *) &ProbconsType, NULL);
}

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC
initprobconslib(void)
{
    PyObject* m;
    
    if (PyType_Ready(&ProbconsType) < 0)
        return;
    
    m = Py_InitModule3("probconslib", NULL, "ProbCons probabilistic consistency based multiple sequence aligner");

    Py_Initialize();
    PyEval_InitThreads();

    Py_INCREF(&ProbconsType);
    PyModule_AddObject(m, "probcons", (PyObject *)&ProbconsType);
    
    char errname[]="Probcons.error";
    ProbconsError = PyErr_NewException(errname, NULL, NULL);
    Py_INCREF(ProbconsError);
    PyModule_AddObject(m, "error", ProbconsError);
}