#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>

const double kB = 1.0;
gsl_rng* rng = NULL;


unsigned long int nFaktorial(unsigned long int n) {
  unsigned long int k, v;
  v = 1;
  for (k = 2; k <= n; k++)
    v *= k;
  return v;
}

unsigned long int nPangkatm(unsigned long int n, unsigned long int m) {
  unsigned long int k, v;
  v = 1;
  for (k = 1; k <= m; k++)
    v *= n;
  return v;
}

int fungsiDelta(long int i, long int j) {
  if (i == j) return 1;
  else return 0;
}


class BitTableField {
/*
  BitTableField class:
  This class provide data structures and generic routines (such
  as e.g. query, insert, delete) for the handling of bit tables.
*/

public:
  unsigned int StaticBitTables = 0;
  unsigned long int static_size = 0;
  int* static_bt = NULL;

  // constructor
  BitTableField (unsigned long int s, bool flag = 0) {

  unsigned long int i;

  IsStatic = flag;

  if (IsStatic) {

    if (static_bt == NULL) {
      static_size = s;
      static_bt = new int[static_size];
      for (i = 0; i < static_size; i++) static_bt[i] = -1;
    }

    StaticBitTables++;
    size = static_size;
    bt = static_bt;

  }

  else {

    size = s;
    bt = new int[size];
    for (i = 0; i < size; i++) bt[i] = -1;

  }

}

  // destructor
  ~BitTableField() {

  if (IsStatic) StaticBitTables--;

  if ((!IsStatic) || (StaticBitTables == 0)) delete[] bt;

  if ((IsStatic) && (StaticBitTables == 0)) {
    static_size = 0;
    static_bt = NULL;
  }

}

  int Query(unsigned long int key)
{

  return bt[key];

}

  void Insert(unsigned long int key, int index)
{

  bt[key] = index;

}
  void Replace(unsigned long int key, int index)
{

  bt[key] = index;

}


  void Delete(unsigned long int key)
{

  bt[key] = -1;

}


  int CheckKey(unsigned long int key)
{

  return bt[key];

}

  unsigned long int GetFillSize()
{

  unsigned long int i, n;

  n = 0;
  for (i = 0; i < size; i++)
    if (bt[i] != -1) n++;

  return n;

}

private:

  int* bt;
  bool IsStatic;
  unsigned long int size;

};

//###### Vector.C ##############################################
//typedef double CoordsType;
typedef long int CoordsType;


// see below...
class Vector;


/*
  Matrix class:

  This class provides a data structure to store m x n matrices of any
  size (dynamically allocated) and various routines for the handling
  of matrix-matrix and matrix-vector operations, e.g. matrix-matrix
  addition/multiplication, matrix-vector multiplication, etc.
*/

class Matrix {

public:

  // constructors:

  // generates a new matrix
  // 1. number of rows
  // 2. number of columns
  // 3. initial value for all matrix elements
  Matrix(int r = 1, int c = 1, CoordsType s = 0)
{
  int i, j;
  rows = r;
  cols = c;
  e = new CoordsType*[rows];
  for (i = 0; i < rows; i++) {
    e[i] = new CoordsType[cols];
    for (j = 0; j < cols; j++)
      e[i][j] = s;
  }
}
  // copy constructor
  Matrix(const Matrix& m)
{
  int i, j;
  rows = m.rows;
  cols = m.cols;
  e = new CoordsType*[rows];
  for (i = 0; i < rows; i++) {
    e[i] = new CoordsType[cols];
    for (j = 0; j < cols; j++)
      e[i][j] = m.e[i][j];
  }
}
  // destructor
  ~Matrix()
{
  int i;
  for (i = 0; i < rows; i++)
    delete[] e[i];
  delete[] e;
}

  // returns a reference to the [i,j]'th matrix element
  // --> "inline" for efficient reading and writing
  CoordsType& Elem(int i, int j) {
    return e[i][j];
  }

  // the subsequent routines are mostly self-explanatory...

  void Set(int i, int j, const CoordsType s)
{
  e[i][j] = s;
}
  
  // special matrix-vector multiplication that performs
  // the "pivot operation"
  // 1. reference to vector which defines the "pivot center"
  // 2. reference to vector on which the operation is performed
  // 3. reference to vector which stores the result
  void PivotOperation(const Vector&, const Vector&, Vector&);

private:

  int rows;        // number of rows
  int cols;        // number of columns
  CoordsType** e;  // pointer to pointer to array with matrix elements

};


/*
  Vector class:

  This class provides a data structure to store vectors of any size
  (dynamically allocated) and various routines for the handling of
  vector operations, e.g. vector-vector addition, scalar product,
  distance between two vectors, etc.
*/

class Vector {

  // These routines are defined in the class "Matrix". They are
  // declared "friend" of the class "Vector" for efficient
  // computation.
  friend void Matrix::PivotOperation(const Vector&, const Vector&, Vector&);

public:

  // constructors:

  // generates a new vector
  // 1. number of vector components
  // 2. initial value for all coordinates
  Vector(int d = 1, CoordsType s = 0)
{
  int i;
  dim = d;
  c = new CoordsType[dim];
  for (i = 0; i < dim; i++)
    c[i] = s;
}
  // copy constructor
  Vector(const Vector& v)
{
  int i;
  dim = v.dim;
  c = new CoordsType[dim];
  for (i = 0; i < dim; i++)
    c[i] = v.c[i];
}
  // destructor
  ~Vector()
{
  delete[] c;
}

  // returns a reference to the i'th vector component
  // --> "inline" for efficient reading and writing
  CoordsType& Elem(int i) {
    return c[i];
  }

  // the subsequent routines are mostly self-explanatory

  CoordsType Get(int i)
{
  return c[i];
}

  void Print()
{
  int i;
  for (i = 0; i < dim; i++)
    printf(" %5ld", c[i]);
  printf("\n");
}

  void Copy(const Vector& v)
{
  int i;
  for (i = 0; i < dim; i++)
    c[i] = v.c[i];
}

  void Copy(const Vector& v1, const Vector& v2)
{
  int i;
  for (i = 0; i < dim; i++)
    c[i] = v2.c[i] - v1.c[i];
}

  CoordsType Distance2(const Vector& v)
{
  int i;
  CoordsType d;
  CoordsType s = 0;
  for (i = 0; i < dim; i++) {
    d = v.c[i] - c[i];
    s += d * d;
  }
  return s;
}

  // checks the overlap between two vectors (on-lattice)
  // returns 1 if overlap, 0 otherwise
  int Overlap(const Vector& v)
{
  int i;
  for (i = 0; i < dim; i++)
    if (v.c[i] != c[i]) return 0;
  return 1;
}

  // checks the overlap between two vectors
  // off-lattice --> 2. threshold
  // returns 1 if overlap, 0 otherwise
  int Overlap(const Vector& v, const CoordsType r)
{
  int i;
  CoordsType d;
  CoordsType s = 0;
  for (i = 0; i < dim; i++) {
    d = v.c[i] - c[i];
    s += d * d;
  }
  if (s <= (r*r)) return 1;
  else return 0;
}

  // returns dimension of first non-equal component between two
  // vectors
  int GetNormalDim(const Vector& v)
{
  int i;
  for (i = 0; i < dim; i++)
    if (v.c[i] != c[i]) return i;
  return -1;
}

private:

  int dim;        // number of vector components
  CoordsType* c;  // pointer to 1D array with vector coordinates

};

void Matrix::PivotOperation(const Vector& p, const Vector& v, Vector& vr)
{
  int i, j;
  for (i = 0; i < vr.dim; i++) {
    vr.c[i] = 0;
    for (j = 0; j < cols; j++)
      vr.c[i] += e[i][j] * (v.c[j] - p.c[j]);
    vr.c[i] += p.c[i];
  }
}


// ###### Histogram.C ##############################################

//typedef double ObservableType;
typedef long int ObservableType;


/*
  Histogram class:

  This class provides data structures and routines for the storing and
  modifying of multi-dimensional histogram data (i.e. DOS, histogram,
  and mask) used for Wang-Landau sampling. Internally, this data is
  stored in dynamically allocated one-dimensional arrays.
*/

class Histogram {

public:

  // constructors:

  // generates a "Histogram" class instance
  // 1. number of dimensions
  // 2. number of physical dimensions
  // 3. array with the boundary specifications for each dimension
  //    i.e. minimum, maximum of the interval and bin size
  Histogram()
{

  int i, j;
  long int n;

  // initialize 'Histogram' class parameters
  Dim = 1;
  DimPhys = 1;
  Bins = 1;
  BinsPhys = 1;
  BinsNonPhys = 1;
  BinsMasked = 0;
  BinsMaskedPhys = 0;
  BinsVisited = 0;
  BinsVisitedPhys = 0;
  Hits = 0;
  HitsPhys = 0;
  YoYo = 1;


  HMin  = new ObservableType[Dim];
  HMax  = new ObservableType[Dim];
  HSize = new ObservableType[Dim];
  HBins = new long int[Dim];
  HInt  = new long int[Dim];

  HMin[0] = 0;
  HMax[0] = 200;
  HSize[0] = 1;

  for (i = 0; i < Dim; i++) {
    HBins[i] = (long int)((HMax[i] - HMin[i]) / HSize[i]);
    // required to keep within array boundaries:
    HMax[i] = HMin[i] + (ObservableType)(HBins[i]) * HSize[i];
    Bins *= HBins[i];
    if (i < DimPhys) BinsPhys *= HBins[i];
  }

  BinsNonPhys = Bins / BinsPhys;

  for (i = Dim - 1; i >= 0; i--) {
    HInt[i] = 1;
    for (j = i + 1; j < Dim; j++)
      HInt[i] *= HBins[j];
  }

  dos  = new double[Bins];
  mask = new int[Bins];
  hist = new unsigned long int[Bins];

  for (n = 0; n < Bins; n++) {
    dos[n] = 0.0;
    mask[n] = 0;
    hist[n] = 0;
  }

}

  // destructor
  ~Histogram()
{

  delete[] dos;
  delete[] mask;
  delete[] hist;

  delete[] HMin;
  delete[] HMax;
  delete[] HSize;
  delete[] HBins;
  delete[] HInt;

}

  // writes the entire histogram state
  // (i.e. specs. and data) to a file
  // 1. a counter
  // 2. file name
  // 3. flag: true = shift DOS to min. DOS; false = don't
  void SaveState(const unsigned long int c = 0, const char* s = NULL, const bool ShiftDOS = true)
{

  int i;
  long int n;
  char filename[50];
  double dosmin;

  if (ShiftDOS) dosmin = GetDOSMin();
  else dosmin = 0.0;

  if (c == 0)
    if (s == NULL) sprintf(filename, "hdata_current.dat");
    else sprintf(filename, "%s", s);
  else
    if (s == NULL) sprintf(filename, "hdata%08lu.dat", c);
    else sprintf(filename, "%s%08lu.dat", s, c);

  FILE* f = fopen(filename, "w");

  fprintf(f, "# State file for 'Histogram' class instance\n");
  fprintf(f, "# Do not edit this file!\n");

  fprintf(f, "Dim  %d\n", Dim);
  fprintf(f, "DimPhys  %d\n", DimPhys);

  fprintf(f, "Bins  %ld\n", Bins);
  fprintf(f, "BinsPhys  %ld\n", BinsPhys);
  fprintf(f, "BinsMasked  %ld\n", BinsMasked);
  fprintf(f, "BinsMaskedPhys  %ld\n", BinsMaskedPhys);
  fprintf(f, "BinsVisited  %ld\n", BinsVisited);
  fprintf(f, "BinsVisitedPhys  %ld\n", BinsVisitedPhys);

  fprintf(f, "Hits  %lu\n", Hits);
  fprintf(f, "HitsPhys  %lu\n", HitsPhys);

  fprintf(f, "YoYo  %d\n", YoYo);

  fprintf(f, "\n");
  fprintf(f, "# Sequence: Dim, HMin, HSize, HBins, HInt\n");
  for (i = 0; i < Dim; i++)
    fprintf(f, " %2d %8ld %8ld %5ld %8ld\n", i, HMin[i], HSize[i], HBins[i], HInt[i]);
  //  fprintf(f, " %2d %12.5e %12.5e %5ld %8ld\n", i, HMin[i], HSize[i], HBins[i], HInt[i]);

  fprintf(f, "\n");
  fprintf(f, "# Sequence: Mask, Density of States (DOS), Histogram\n");
  for (n = 0; n < Bins; n++)
    if (mask[n])
      fprintf(f, " %d  %.16e  %lu\n", mask[n], dos[n] - dosmin, hist[n]);
    else
      fprintf(f, " %d  %.16e  %lu\n", mask[n], dos[n], hist[n]);

  fclose(f);

}

  // returns min. and max. of the DOS
  // 1. flag: 0 = all dim.; 1 = only physical dim.
  double GetDOSMin(int OnlyPhys = 0)
{

  long int n, incr;
  int flag = 1;
  double min = 0.0;

  if (OnlyPhys) incr = BinsNonPhys;
  else incr = 1;

  for (n = 0; n < Bins; n += incr)
    if (mask[n]) {
      if (flag) {
	flag = 0;
	min = dos[n];
      }
      else if (dos[n] < min) min = dos[n];
    }

  return min;

}
  
  double GetDOSMax(int OnlyPhys = 0)
{

  long int n, incr;
  int flag = 1;
  double max = 0.0;

  if (OnlyPhys) incr = BinsNonPhys;
  else incr = 1;

  for (n = 0; n < Bins; n += incr)
    if (mask[n]) {
      if (flag) {
	flag = 0;
	max = dos[n];
      }
      else if (dos[n] > max) max = dos[n];
    }

  return max;

}

  // prints normalized (to unity) DOS
  // 1. file name
  // 2. flag: 0 = all dim.; 1 = only physical dim.
  void PrintNormDOS(const char* filename = NULL, int OnlyPhys = 1)
{

  const double dosmax = GetDOSMax(OnlyPhys);

  long int m, n, incr;
  int d, i;
  double c = 0.0;

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  if (OnlyPhys) {
    incr = BinsNonPhys;
    d = DimPhys;
  }
  else {
    incr = 1;
    d = Dim;
  }

  for (n = 0; n < Bins; n += incr)
    if (mask[n]) c += exp(dos[n] - dosmax);

  for (n = 0; n < Bins; n += incr)
    if (mask[n]) {
      m = n;
      for (i = 0; i < d; i++) {
	fprintf(f, " %4ld", HMin[i] + (m / HInt[i]) * HSize[i]);
	m %= HInt[i];
      }
      fprintf(f, "  %18.10e\n", exp(dos[n] - dosmax) / c);
    }

  if (filename != NULL) fclose(f);

}

  // checks the histogram for Wang-Landau sampling
  // 1. check type: 0 = flatness check (flatness criterion)
  //                1 = minimal number of hits per bin
  // 2. flatness criterion (0 - 1)
  // 3. minimal number of hits per bin
  int CheckHistogram(int CheckType, double p, unsigned long int bmin)
{

  long int n;
  double threshold;

  switch (CheckType) {

  case 0 : {

    threshold = p * (double)(Hits) / (double)(BinsMasked);

    for (n = 0; n < Bins; n++)
      if ((mask[n]) && ((double)(hist[n]) < threshold)) return 0;

    return 1;

  }

  case 1 : {

    for (n = 0; n < Bins; n++)
      if ((mask[n]) && (hist[n] < bmin)) return 0;

    return 1;

  }

  }

  return 0;

}


  // resets the histogram
  void ResetHist()
{
  long int n;
  for (n = 0; n < Bins; n++)
    hist[n] = 0;
  BinsVisited = 0;
  BinsVisitedPhys = 0;
  Hits = 0;
  HitsPhys = 0;
}

  // returns the histogram array index or -1 if out of boundaries
  // 1. 1D array of observables
  // 2. flag: 0 = all bins; 1 = only within masked bins
  long int GetIndex(ObservableType v[], int MaskCheck = 0)
{

  int i;
  long int index = 0;

  for (i = 0; i < Dim; i++) {
    if ((v[i] < HMin[i]) || (v[i] >= HMax[i])) return -1;
    index += HInt[i] * (long int)((v[i] - HMin[i]) / HSize[i]);
  }

  if ((MaskCheck) && (!mask[index])) return -1;
  else return index;

}

  // updates the DOS, histogram and mask for Wang-Landau sampling
  // 1. histogram array index
  // 2. reference to the counter of the histogram check interval
  // 3. update type
  // 4. DOS modification factor
  void Update(long int index,
		       unsigned long int& counter,
		       int UpdateType, double f)
{

  int i, PhysFlag;

  if ((index % BinsNonPhys) == 0) PhysFlag = 1;
  else PhysFlag = 0;

  switch (UpdateType) {

  case  1 : {
    if (mask[index]) dos[index] += f;
    else {
      ResetHist();
      counter = 0;
      dos[index] = GetDOSMin();
    }
    break;
  }

  default : {
    dos[index] += f;
  }

  }

  if (!hist[index]) {
    BinsVisited++;
    if (PhysFlag) BinsVisitedPhys++;
  }

  hist[index]++;
  Hits++;
  if (PhysFlag) HitsPhys++;

  if (!mask[index]) {
    mask[index] = 1;
    BinsMasked++;
    if (PhysFlag) BinsMaskedPhys++;
    printf("New bin masked:");
    for (i = 0; i < Dim; i++) {
      printf(" %4ld", HMin[i] + (index / HInt[i]) * HSize[i]);
      index %= HInt[i];
    }
    printf("\n");
    printf("Actual DOS minimum  = %15.8e\n", GetDOSMin());
    printf("Actual mod. factor  = %15.8e\n", f);
    printf("Ratio (Min DOS / f) = %15.8e\n", GetDOSMin() / f);
  }

}

  // gets DOS, histogram and mask ("inline")
  // 1. histogram array index
  double GetDOS(const long int index) {
    return dos[index];
  }
  int GetMask(const long int index) {
    return mask[index];
  }
  unsigned long int GetHist(const long int index) {
    return hist[index];
  }

  // records and prints the "up" and "down" round trip times
  // returns 1 in case of an up-walk, otherwise returns 0
  // 1. histogram array index
  // 2. current number of Monte Carlo steps
  // 3. reference to a variable that stores the MC steps
  // 4. dimension for which the round trips are to be checked
  // 5. file name (optional)
  int CheckItinerancy(long int index,
			       const unsigned long int MCMoves,
			       unsigned long int& MCMovesMem,
			       int dim = 0, const char* filename = NULL)
{

  int i;
  long int n;
  FILE* f;

  if ((index % BinsNonPhys) == 0) {  // only "physical" round-trips are considered

    for (i = 0; i < dim; i++) index %= HInt[i];
    index /= HInt[dim];

    if (YoYo == 1) {  // up-walk

      for (n = Bins - BinsNonPhys; n >= 0; n -= BinsNonPhys)
	if (mask[n]) {
	  for (i = 0; i < dim; i++) n %= HInt[i];
	  n /= HInt[dim];
	  break;
	}

      if (index == n) {
	if (filename != NULL) f = fopen(filename, "a");
	else f = stdout;
	fflush(f);
	if (filename != NULL) fclose(f);
	MCMovesMem = MCMoves;
	YoYo = 0;
	return 1;
      }

    }

    else {  // down-walk

      for (n = 0; n < Bins; n += BinsNonPhys)
	if (mask[n]) {
	  for (i = 0; i < dim; i++) n %= HInt[i];
	  n /= HInt[dim];
	  break;
	}

      if (index == n) {
	if (filename != NULL) f = fopen(filename, "a");
	else f = stdout;
	fflush(f);
	if (filename != NULL) fclose(f);
	MCMovesMem = MCMoves;
	YoYo = 1;
      }

    }

  }

  return 0;

}

private:

  int Dim;      // total number of dimensions
  int DimPhys;  // number of physical dimensions

  long int Bins;             // total number of bins
  long int BinsPhys;         // number of physical bins
  long int BinsNonPhys;      // Bins / BinsPhys
  long int BinsMasked;       // total number of masked bins
  long int BinsMaskedPhys;   // number of physical masked bins
  long int BinsVisited;      // total number of currently visited bins
  long int BinsVisitedPhys;  // number of currently physical visited bins

  unsigned long int Hits;      // sum of histogram entries
  unsigned long int HitsPhys;  // sum of physical histogram entries

  ObservableType* HMin;   // pointer to 1D array with histogram boundary minima
  ObservableType* HMax;   // pointer to 1D array with histogram boundary maxima
  ObservableType* HSize;  // pointer to 1D array with histogram bin sizes

  long int* HBins;  // pointer to 1D array with histogram sizes (for each dim.)
  long int* HInt;   // pointer to 1D array with histogram "size offsets"

  int* mask;                // pointer to the 'mask' (1D array)
  double* dos;              // pointer to the DOS (1D array)
  unsigned long int* hist;  // pointer to the histogram (1D array)

  int YoYo;  // flag used by CheckItinerancy: 0 = "down walk"; 1 = "up walk"

};



// ###### H0PModel.C ##############################################
/*
  H0PModel class:

  This is a  class that contains all relevant data structures and routines to
  run Monte Carlo simulations for the H0P protein
  model on the simple (hyper-)cubic lattice in dimension d >= 2. For efficiency
  reasons (e.g. performing Monte Carlo trial moves, checking
  self-overlaps, counting non-bonded HH contacts), the configuration
  of the protein is stored both as a set of d-dimensional
  vectors for the monomer coordinates and as entries in a
  d-dimensional occupancy field. Currently, the following three types
  of Monte Carlo trial moves are implemented: pull moves,
  bond-rebridging moves, pivot moves. The principal observable is the
  number of non-bonded HH contacts (the energy of these models is just
  the negative of this quantity). To each Monte Carlo trial move
  belongs the updating of this observable.
*/
class H0PModel {

public:

  // constructor
  // generates a "H0PModel" class instance
  // 1. input file name
  H0PModel(const char* filename)
{

  const int MaxSeqH0P = 200000;
  char line[200];
  FILE* f;
  
  char s[200];
  char* sH0P;

  int i, n;
  char c;

  f = fopen(filename, "r");
  fgets(line, sizeof(line), f);
  sscanf(line, "%s", s);

  
  // initialize 'H0PModel' class parameters
  PolymerDim = 3;
  OccupancyFieldType = 0;
  MoveFraction[0] = 0.75;
  MoveFraction[1] = 0.98;

  NumberOfMonomers = 0;
  sH0P = new char[MaxSeqH0P + 1];
  sH0P[0] = '\0';

	if ((int)(strlen(sH0P) + strlen(s)) <= MaxSeqH0P) {
	  strcpy(sH0P + strlen(sH0P), s);
	  NumberOfMonomers = strlen(sH0P);
  }
  
  // now, construct class instance...

  // pull-2 moves may shift the polymer by max. 2 lattice units
  // --> LSize must be >= polymer chain length + 2 in order to
  //     allow for a "free moving" of the polymer in space
  LSize = NumberOfMonomers + 2;
  LOffset = new unsigned long int[PolymerDim];
  for (n = 0; n < PolymerDim; n++)
    LOffset[n] = nPangkatm(LSize, PolymerDim - 1 - n);

  seq = new MonomerType[NumberOfMonomers];
  coords = new Vector*[NumberOfMonomers];
  tmpcoords = new Vector*[NumberOfMonomers];

  for (n = 0; n < NumberOfMonomers; n++) {
    coords[n] = new Vector(PolymerDim);
    tmpcoords[n] = new Vector(PolymerDim);
    coords[n] -> Elem(0) = n;
    switch (sH0P[n]) {
    case 'P' : {seq[n] = polar; break;}
    case 'H' : {seq[n] = hydrophobic; break;}
    case '0' : {seq[n] = neutral; break;}
    default  : {}
    }
  }

  // read coordinates file...
  if (f == NULL) f = fopen("coords_init.xyz", "r");

  if (OccupancyFieldType == 0)
    OccupancyField = new BitTableField(nPangkatm(LSize, PolymerDim));

  site = new unsigned long int[NumberOfMonomers];
  tmpsite = new unsigned long int[NumberOfMonomers];

  for (n = 0; n < NumberOfMonomers; n++) {
    site[n] = GetSite(*coords[n]);
    tmpsite[n] = site[n];
    OccupancyField -> Insert(site[n], n);
  }

  NumberOfHMonomers = 0;
  for (n = 0; n < NumberOfMonomers; n++)
    if (seq[n] == hydrophobic) NumberOfHMonomers++;

  InternalHHContacts = 0;
  for (n = 0; n < NumberOfMonomers - 1; n++)
    if ((seq[n] == hydrophobic) &&
	(seq[n + 1] == hydrophobic)) InternalHHContacts++;

  InternalH0Contacts = 0;
  for (n = 0; n < NumberOfMonomers - 1; n++) {
    if ((seq[n] == hydrophobic) && (seq[n + 1] == neutral)) InternalH0Contacts++;
    if ((seq[n] == neutral) && (seq[n + 1] == hydrophobic)) InternalH0Contacts++;
  }

  CheckModelIntegrity();

  // pull-2 move:
  InitPull2Set();

  // bond-rebridging move: for relabeling/restoring
  CoordsPointer = new Vector*[NumberOfMonomers];
  CoordsIndex = new int[NumberOfMonomers];

  // pivot move:
  InitPivotMatrices();

  // observables:
  // measured during Monte Carlo sampling
  // 0  : number of non-bonded HH contacts (principal observable)

  // additional observables could be defined and sampled, e.g. for a
  // 2D Wang-Landau simulation

  InitObservables(1);

  // evaluate "initial" principal observables...

  // H0P model Hamiltonian
  H0PModelFlag = 1;
  GetObservablesH0P();
  
  // --> if (BeadsInterval < BeadsIntervalThreshold) use 'GetObservablesPart'
  //     otherwise use 'GetObservables'
  BeadsIntervalThreshold = NumberOfMonomers / 4;

  // required for Monte Carlo moves statistics
  for (i = 0; i < 3; i++) {
    MoveTrials[i] = 0;
    MoveRejects[i] = 0;
    MoveDiscards[i] = 0;
    BeadsTrials[i] = 0;
    BeadsRejects[i] = 0;
  }

  WriteState(0);

  delete[] sH0P;

}

  // destructor
  ~H0PModel()
{

  int n;
  long int i;

  CheckModelIntegrity();

  delete[] seq;

  for (n = 0; n < NumberOfMonomers; n++) {
    delete coords[n];
    delete tmpcoords[n];
  }
  delete[] coords;
  delete[] tmpcoords;

  delete OccupancyField;

  delete[] site;
  delete[] tmpsite;

  delete[] Pull2Set;

  delete[] CoordsPointer;
  delete[] CoordsIndex;


  DeleteObservables();

  for (i = 0; i < NPiOps; i++)
    delete pivot[i];
  delete[] pivot;

  delete[] LOffset;

  WriteState(5);

  printf("Ciao H0PModel...\n");

}

  // see comments for the basis class "Model" for the following three
  // routines
  void WriteState(int format = 0, const char* filename = NULL)
{

  const char symbol[3] = {'H', 'P' , '0'};

  int i, n;

  FILE* f;
  if (filename != NULL) f = fopen(filename, "a");
  else f = stdout;

  switch (format) {

    // mol2 format output
  case  1 : {

    double* ctr = new double[PolymerDim];
    for (i = 0; i < PolymerDim; i++) ctr[i] = 0.0;
    for (n = 0; n < NumberOfMonomers; n++)
      for (i = 0; i < PolymerDim; i++)
	ctr[i] += (double)(coords[n] -> Elem(i));
    for (i = 0; i < PolymerDim; i++) ctr[i] /= (double)(NumberOfMonomers);

    fprintf(f, "\n");
    fprintf(f, "@<TRIPOS>MOLECULE\n");
    fprintf(f, "H0PMODEL\n");
    fprintf(f, " %d %d\n", NumberOfMonomers, NumberOfMonomers - 1);
    fprintf(f, "PROTEIN\n");
    fprintf(f, "NO_CHARGES\n");
    fprintf(f, "@<TRIPOS>ATOM\n");
    for (n = 0; n < NumberOfMonomers; n++) {
      fprintf(f, "%5d %c", n + 1, symbol[seq[n]]);
      if (PolymerDim == 2) {
	fprintf(f, " %9.3f %9.3f %9.3f",
		(double)(coords[n] -> Elem(0)) - ctr[0],
		(double)(coords[n] -> Elem(1)) - ctr[1],
		0.0);
      }
      else {
	fprintf(f, " %9.3f %9.3f %9.3f",
		(double)(coords[n] -> Elem(0)) - ctr[0],
		(double)(coords[n] -> Elem(1)) - ctr[1],
		(double)(coords[n] -> Elem(2)) - ctr[2]);
      }
      fprintf(f, " %c\n", symbol[seq[n]]);
    }
    fprintf(f, "@<TRIPOS>BOND\n");
    for (n = 1; n < NumberOfMonomers; n++)
      fprintf(f, "%5d %5d %5d 1\n", n, n, n + 1);

    delete[] ctr;

    break;

  }

    // pdb format output
  case  2 : {

    fprintf(f, "MODEL\n");
    for (n = 0; n < NumberOfMonomers; n++) {
      if (seq[n] == polar)
	fprintf(f, "ATOM   %4d  O           0", n + 1);
      else
	fprintf(f, "ATOM   %4d  P           0", n + 1);
      if (PolymerDim == 2) {
	fprintf(f, " %9.2f %9.2f %9.2f\n",
		10.0 * (double)(coords[n] -> Get(0) - coords[0] -> Get(0)),
		10.0 * (double)(coords[n] -> Get(1) - coords[0] -> Get(1)),
		0.0);
      }
      else {
	fprintf(f, " %9.2f %9.2f %9.2f\n",
		10.0 * (double)(coords[n] -> Get(0) - coords[0] -> Get(0)),
		10.0 * (double)(coords[n] -> Get(1) - coords[0] -> Get(1)),
		10.0 * (double)(coords[n] -> Get(2) - coords[0] -> Get(2)));
      }
    }
    fprintf(f, "CONECT %4d %4d\n", 1, 2);
    for (n = 2; n < NumberOfMonomers; n++)
      fprintf(f, "CONECT %4d %4d %4d\n", n, n - 1, n + 1);
    fprintf(f, "CONECT %4d %4d\n", NumberOfMonomers, NumberOfMonomers - 1);
    fprintf(f, "ENDMDL\n\n");

    break;

  }

    // current model state and observables output
  case  3 : {

    fprintf(f, "\n");
    fprintf(f, "H0P/ISAW Model : %d\n", NumberOfMonomers);
    PrintObservables(f);
    for (n = 0; n < NumberOfMonomers; n++) {
      fprintf(f, "%c", symbol[seq[n]]);
      for (i = 0; i < PolymerDim; i++)
	fprintf(f, " %ld", coords[n] -> Get(i) - coords[0] -> Get(i));
      fprintf(f, "\n");
    }
    break;

  }

    // write current model state
  case  4 : {

    if (filename != NULL) {
      fclose(f);
      f = fopen(filename, "w");
    }
    for (n = 0; n < NumberOfMonomers; n++) {
      fprintf(f, "%c", symbol[seq[n]]);
      for (i = 0; i < PolymerDim; i++)
	fprintf(f, " %ld", coords[n] -> Get(i));
      fprintf(f, "\n");
    }
    break;

  }

    // print move statistics
  case  5 : {

    fprintf(f, "\nMove statistics:\n\n");
    fprintf(f, " Move type       | # Attempts   | Valid moves | Acceptance ratio  \n");
    fprintf(f, "                 |              | (ratio)     | (valid moves)     \n");
    fprintf(f, "-----------------|--------------|-------------|------------------ \n");
    for (i = 0; i < 3; i++) {
      switch (i) {
      case 0 : {fprintf(f, " pull-2          |"); break;}
      case 1 : {fprintf(f, " bond-rebridging |"); break;}
      case 2 : {fprintf(f, " pivot           |"); break;}
      }
      fprintf(f, " %12ld |  %10.4e |       %10.4e\n",
	      MoveTrials[i] + MoveDiscards[i],
	      double(MoveTrials[i]) / double(MoveTrials[i] + MoveDiscards[i]),
	      double(MoveTrials[i] - MoveRejects[i]) / double(MoveTrials[i]));
    }
    fprintf(f, "\n");
    break;

  }

    // print model specifications
  case  0 : {

    printf("\n### Physical Model #####################################\n\n");
    printf("H0P Model in d dimensions\n\n");
    printf("Specifications:\n\n");
    printf(" Dimensions = %d\n", PolymerDim);
    printf(" Number of monomers = %d (%d)\n", NumberOfMonomers, NumberOfHMonomers);
    printf(" Internal HH contacts = %d\n", InternalHHContacts);
    printf(" Observables type = ");
    if (H0PModelFlag) printf("H0P model\n");
    else printf("ISAW model\n");
    printf(" Polymer fragment threshold = %d\n", BeadsIntervalThreshold);
    printf(" Lattice size = %ld (%lu)\n", LSize, nPangkatm(LSize, PolymerDim));
    printf(" Move fractions = %6.4f / %6.4f / %6.4f\n",
	   MoveFraction[0],
	   MoveFraction[1] - MoveFraction[0],
	   1.0 - MoveFraction[1]);
    printf("\nInitial state:\n\n");
    for (n = 0; n < NumberOfMonomers; n++) {
      printf(" %c ", symbol[seq[n]]);
      coords[n] -> Print();
    }
    printf("\n########################################################\n");
    break;

  }

  }

  if (filename != NULL) fclose(f);

}

  void DoMCMove()
{

  double r = gsl_rng_uniform(rng);

  // backup "old" principal observable...
  tmp_observable[0] = observable[0];

  MoveProposal = 0;

  if (r < MoveFraction[0]) {       // pull-2 move
    MCMoveType = 0;
    Pull2Move();
  }

  else if (r < MoveFraction[1]) {  // bond-rebridging move
    MCMoveType = 1;
    BondRebridgingMove();
  }

  else {                           // pivot move
    MCMoveType = 2;
    PivotMove();
  }

  if (MoveProposal) {
    MoveTrials[MCMoveType]++;
    BeadsTrials[MCMoveType] += BeadsInterval;
  }
  else MoveDiscards[MCMoveType]++;

  //  printf("Do MC move %d ...\n", MCMoveType);
  //  CheckModelIntegrity();

}

  void UnDoMCMove()
{

  int n;
  Vector* t;
  Vector** tt;
  unsigned long int* s;

  //  printf("Undo MC move %d ...\n", MCMoveType);

  MoveRejects[MCMoveType]++;
  BeadsRejects[MCMoveType] += BeadsInterval;

  // restore "old" principal observable...
  observable[0] = tmp_observable[0];

  // pull-2 move
  if (MCMoveType == 0) {

    n = BeadsEnd;

    while (1) {

      OccupancyField -> Delete(site[n]);
      site[n] = tmpsite[n];
      OccupancyField -> Insert(site[n], n);

      t = coords[n];
      coords[n] = tmpcoords[n];
      tmpcoords[n] = t;

      if (n == BeadsStart) break;

      n -= BeadsIncr;

    }

  }

  // bond-rebridging move
  else if (MCMoveType == 1) {

    for (n = 0; n < NumberOfMonomers; n++)
      if (CoordsIndex[n] != n)
	OccupancyField -> Replace(site[n], CoordsIndex[n]);

    tt = coords;
    coords = CoordsPointer;
    CoordsPointer = tt;

    s = site;
    site = tmpsite;
    tmpsite = s;

  }

  // pivot move
  else {

    n = BeadsStart;
    while (1) {
      OccupancyField -> Delete(site[n]);
      t = coords[n];
      coords[n] = tmpcoords[n];
      tmpcoords[n] = t;
      if (n == BeadsEnd) break;
      n += BeadsIncr;
    }

    n = BeadsStart;
    while (1) {
      site[n] = tmpsite[n];
      OccupancyField -> Insert(site[n], n);
      if (n == BeadsEnd) break;
      n += BeadsIncr;
    }

  }

  //  CheckModelIntegrity();

}


  // number of defined observables for the "physical model"; for the
  // "H0PModel" class, there is currently only one observable defined
  // (number of non-bonded HH contacts)
  int NumberOfObservables;
  // pointer to the array of observables
  ObservableType* observable;

  // state variable indicating whether a particular Monte Carlo trial
  // move has been applied successfully (i.e. no self-overlap
  // violation), (0 = not successful; 1 = successful); in case of not
  // successful, the time consuming Wang-Landau / Metropolis
  // acceptance criterion check can be skipped
  int MoveProposal;

  ObservableType* tmp_observable;

  void PrintObservables(FILE* f = stdout) {

    int i;

    fprintf(f, "Observables:");
    for (i = 0; i < NumberOfObservables; i++)
      fprintf(f, " %ld", observable[i]);
    fprintf(f, "\n");

  }

  void InitObservables(int n) {

    int i;

    NumberOfObservables = n;

    observable = new ObservableType[NumberOfObservables];
    tmp_observable = new ObservableType[NumberOfObservables];

    for (i = 0; i < NumberOfObservables; i++) {
      observable[i] = (ObservableType)(0);
      tmp_observable[i] = (ObservableType)(0);
    }

  }

  void DeleteObservables() {

    delete[] observable;
    delete[] tmp_observable;

  }
  
  double dummyRg[100];
  double Rg[100];
  long int Rgcount[100];
  double COM[3];

  void GetRadiusofGyration() {
    // initiating variable to calculate radius of gyration
    
    double currRg;

    
    // calculating center of mass
    for (int i = 0; i < 3; i++) {
      for (int n = 0; n < NumberOfMonomers; n++) {
        double tmpCOM = coords[n] -> Elem(i);
        COM[i] += tmpCOM/NumberOfMonomers;        
      }
    }
    /*
    // print out the center of mass of the current structure
    for (int i = 0; i < 3; i++) {
    printf("COM = %f \n", COM[i]);
    fflush(stdout);
    }
    */

    // calculating radius of gyration for current structure
    double gx = 0, gy = 0, gz = 0;

    for (int n = 0; n < NumberOfMonomers; n++) {
      gx += (coords[n] -> Elem(0)  - COM[0])*(coords[n] -> Elem(0)  - COM[0]);
      gy += (coords[n] -> Elem(1)  - COM[1])*(coords[n] -> Elem(1)  - COM[1]);
      gz += (coords[n] -> Elem(2)  - COM[2])*(coords[n] -> Elem(2)  - COM[2]);
    }
    /*  
    printf("gx= %f \n", gx);
    printf("gy= %f \n", gy);
    printf("gz= %f \n", gz);
    fflush(stdout);
    */
    double h = gx + gy + gz;
    currRg = sqrt(h/NumberOfMonomers);
    // printf("currRg = %f \n", currRg);

    // adding the calculated Rg to sum up at the respective energy level
    // printf("nHH = %ld \n", observable[0]);
    dummyRg[observable[0]] += currRg;

    // count the number of structure sampled for each energy level
    Rgcount[observable[0]]++;

    // restarting the center of mass value
    COM[0] = 0.0;
    COM[1] = 0.0;
    COM[2] = 0.0;
    gx = 0;
    gy = 0;
    gz = 0;
      
  }

  void GetAverageRg() {
    // calculate the average radius of gyration for each energy level and print the value
    for (int j = 0; j < 100; j++) {
      Rg[j] = dummyRg[j]/Rgcount[j];
      printf("Rg = %f \n", Rg[j]);
      fflush(stdout);
    }
  }

  int CrossProductMagnitude(int vect_A[], int vect_B[]) {
    int cross_P[3];
    int m;
    cross_P[0] = vect_A[1]*vect_B[2] - vect_A[2]*vect_B[1];
    cross_P[1] = vect_A[2]*vect_B[0] - vect_A[0]*vect_B[2];
    cross_P[2] = vect_A[0]*vect_B[1] - vect_A[1]*vect_B[0];
    m = cross_P[0] + cross_P[1] + cross_P[2];
    return m;
  }
  
  double dummyTau[100];
  double Tau[100];
  long int Taucount[100];

  void GetTortuosity() {
    

    long int Si[NumberOfMonomers-2];
    double Sav;
    double sumSi = 0.0;
    double currTau;
    double SiS = 0.0;
    for (int i = 0; i < NumberOfMonomers-2; i++) {
    long int currSi = 0;
      // calculate the si
      for (int j = 0; j <= i ; j++) {
        int rjj1[3], rjj2[3];
        // calculate the rjj1 and rjj2 vector
        for (int x = 0; x < 3; x++) {
          rjj1[x] = (coords[j+1] -> Elem(x))-(coords[j] -> Elem(x));
          rjj2[x] = (coords[j+2] -> Elem(x))-(coords[j] -> Elem(x));
        }
        currSi += CrossProductMagnitude(rjj1, rjj2);
      }
      Si[i] = currSi;
      sumSi += Si[i];
      currSi = 0;
    }
    Sav = sumSi / (NumberOfMonomers-2);
    for (int i = 0; i < (NumberOfMonomers-2); i++) {
      SiS += (Si[i] - Sav)*(Si[i] - Sav);
    }
    currTau = sqrt(SiS/(NumberOfMonomers-2));
    
    dummyTau[observable[0]] += currTau;
    Taucount[observable[0]]++;

    sumSi = 0.0;
    Sav = 0;
    SiS = 0.0;
    for (int i = 0; i < NumberOfMonomers-2; i++) {
      Si[i] = 0;
    }
  }

  void GetAverageTau() {
    // calculate the average radius of gyration for each energy level and print the value
    for (int j = 0; j < 100; j++) {
      Tau[j] = dummyTau[j]/Taucount[j];
      printf("Tau = %f \n", Tau[j]);
      fflush(stdout);
    }
  }

  void printRg() {
    FILE* f;
    f = fopen("Rg.dat", "w");
    for (int i = 0; i < 100; i++) {
      fprintf(f, " %4d", i);
      fprintf(f, "    %f \n", Rg[i]);
    }
  }

  void printTau() {
    FILE* f;
    f = fopen("Tau.dat", "w");
    for (int i = 0; i < 100; i++) {
      fprintf(f, " %4d", i);
      fprintf(f, "    %f \n", Tau[i]);
    }
  }

private:

  // the two monomer types: "H" and "P"
  enum MonomerType {hydrophobic = 0, polar = 1, neutral = 2};

  // physical dimension d of the protein/polymer (d >= 2)
  int PolymerDim;
  // the occupancy field can be either a bit table or a hash table,
  // see class "Field"
  int OccupancyFieldType;
  // array specifying the relative frequency of each type of Monte
  // Carlo trial move; as an example: [0.5, 0.8] means 50% pull moves,
  // 30% bond-rebridging moves and 20% pivot moves
  double MoveFraction[2];

  // various arrays to collect statistics for the three (3) types of
  // Monte Carlo trial moves
  unsigned long int MoveTrials[3];
  unsigned long int MoveRejects[3];
  unsigned long int MoveDiscards[3];
  unsigned long int BeadsTrials[3];
  unsigned long int BeadsRejects[3];

  // two parameters specifying the d-dimensional occupancy field
  long int LSize;
  unsigned long int* LOffset;

  // total number of pivot operations in d dimensions and pointer to
  // the corresponding pivot matrices
  long int NPiOps;
  Matrix** pivot;

  // total number of monomers (protein/polymer chain length)
  int NumberOfMonomers;
  // number of "H" monomers (equals 'NumberOfMonomers' in case of ISAW
  // model)
  int NumberOfHMonomers;
  // number of chain-internal "HH" contacts (equals
  // 'NumberOfMonomers-1' in case of ISAW model)
  int InternalHHContacts;

  int InternalH0Contacts;

  // flag indicating whether the current system is an "H0P model" or an
  // "ISAW model"; it is an "H0P model" if there is at least one "P"
  // monomer in the sequence; 0 = "ISAW model", 1 = "H0P model"
  int H0PModelFlag;

  // pointer to "Field" class instance (either "BitTable" or
  // "HashTable")
  BitTableField* OccupancyField;

  // pointer to sequence array of "H"'s and "P"'s
  MonomerType* seq;

  // pointer to array of vectors for the monomer coordinates
  Vector** coords;
  // same for temporary storage
  Vector** tmpcoords;

  // pointers to arrays of monomer locations in the d-dimensional
  // occupancy field
  unsigned long int* site;
  unsigned long int* tmpsite;

  // index for the type of Monte Carlo trial move; 0 = pull move; 1 =
  // bond-rebridging move; 2 = pivot move
  int MCMoveType;

  // "first" and "last" monomers (and "direction") of the chain
  // fragment whose monomer positions have changed after a Monte Carlo
  // trial move (only used for pull moves and pivot moves)
  int BeadsStart, BeadsEnd, BeadsIncr;
  // length of chain fragment whose monomer positions have changed
  // after a Monte Carlo trial move (for all types of Monte Carlo
  // trial moves); "threshold" length of chain fragment determining
  // which method to use to evaluate the number of non-bonded HH
  // contacts (either for the entire chain or only for the chain
  // fragment)
  int BeadsInterval, BeadsIntervalThreshold;

  // maximal number of all possible pull moves
  long int MaxPull2Moves;
  // pointer to a "look-up" table with the relative displacements for
  // all possible pull moves
  int (*Pull2Set)[6];

  // pointer to temporary "buffer array" to keep vector addresses for
  // the monomer coordinates
  Vector** CoordsPointer;
  // pointer to temporary "buffer array" to keep monomer indices
  int* CoordsIndex;

  // initializes the pivot matrices in d dimensions
  // in 2D the 8 pivot matrices generated by InitPivotMatrices()
// correspond to the following symmetry operations
// on the square lattice (relative to a given pivot point):

// pivot[0] = identity
// pivot[1] = reflexion on principal diagonal
// pivot[2] = reflexion on y-axis
// pivot[3] = +90 degree rotation (counter-clockwise)
// pivot[4] = reflexion on x-axis
// pivot[5] = -90 degree rotation (clockwise)
// pivot[6] = 180 degree rotation
// pivot[7] = reflexion on secondary diagonal

// in d dimensions, the total number of
// pivot operations is 2^d * d!

void InitPivotMatrices()
{

  const long int c = nPangkatm(2, PolymerDim);

  long int k, n, v;
  int i, j;
  Matrix* t;

  long int* s = new long int[PolymerDim];
  gsl_permutation* p = gsl_permutation_alloc(PolymerDim);

  printf("\n\nInitializing Pivot matrices...");

  NPiOps = c * nFaktorial(PolymerDim);
  pivot = new Matrix*[NPiOps];
  for (n = 0; n < NPiOps; n++)
    pivot[n] = new Matrix(PolymerDim, PolymerDim);

  n = 0;
  for (k = 0; k < c; k++) {

    v = k;
    for (i = PolymerDim - 1; i >= 0; i--) {
      s[i] = v / nPangkatm(2, i);
      s[i] = 1 - 2 * s[i];
      v %= nPangkatm(2, i);
    }

    gsl_permutation_init(p);

    do {

      for (i = 0; i < PolymerDim; i++)
	for (j = 0; j < PolymerDim; j++)
	  pivot[n] ->
	    Set(i, j, s[i] * fungsiDelta(i, gsl_permutation_get(p, j)));

      //      pivot[n] -> Print();
      n++;

    } while (gsl_permutation_next(p) == GSL_SUCCESS);

  }

  if (n != NPiOps) ;

  printf("done\n");
  printf("--> number of Pivot matrices = %ld\n\n", NPiOps);

  gsl_permutation_free(p);
  delete[] s;

  // identity operation (0) does not generate
  // a "new" configuration
  // --> shift to the end of array of pivot matrices
  t = pivot[0];
  pivot[0] = pivot[n - 1];
  pivot[n - 1] = t;

}
  // generates a "look-up" table with the relative displacements for
  // all possible pull moves (labeled "pull-2" moves because there are
  // generally two monomers which are pulled at a time)
  void InitPull2Set()
{

  int n, d, i, u, j, v;
  long int m = 0;

  // max. number of pull-2 moves:
  // --> 2 * (max. number of end pull-2 moves)
  //     + (n - 2) * (max. number of internal pull-2 moves)

  // note: the square lattice in d dimensions has
  //       2 * d nearest neighbor sites and
  //       2 * (d - 1) nearest neighbor sites in the hyper-plane
  //       perpendicular to a given direction

  printf("\n\nInitializing pull-2 moves...");

  MaxPull2Moves = 0;

  // internal pull-2 moves
  MaxPull2Moves += (NumberOfMonomers-2) * 2 * 2*(PolymerDim-1);

  // end pull-2 moves
  MaxPull2Moves += 2 * (2*(PolymerDim-1) + 1);
  MaxPull2Moves += 2 * (2*(PolymerDim-1) * (2*(PolymerDim-2) + 1 + 1));

  Pull2Set = new int[MaxPull2Moves][6];

  for (n = 1; n < NumberOfMonomers - 1; n++) {
    for (d = -1; d <= 1; d += 2) {
      for (i = 1; i <= PolymerDim - 1; i++) {
	for (u = -1; u <= 1; u += 2) {

	  Pull2Set[m][0] = n;
	  Pull2Set[m][1] = d;
	  Pull2Set[m][2] = i;
	  Pull2Set[m][3] = u;
	  Pull2Set[m][4] = -1;
	  Pull2Set[m][5] = -1;

// 	  printf("set %3ld : %3d %3d %3d %3d %3d %3d\n",
// 		 m,
// 		 Pull2Set[m][0],
// 		 Pull2Set[m][1],
// 		 Pull2Set[m][2],
// 		 Pull2Set[m][3],
// 		 Pull2Set[m][4],
// 		 Pull2Set[m][5]);

	  m++;

	}
      }
    }
  }

  for (n = 0; n < 2; n++) {
    for (i = 0; i <= PolymerDim - 1; i++) {
      for (u = -1; u <= 1; u += 2) {
	for (j = 0; j <= PolymerDim - 1; j++) {
	  for (v = -1; v <= 1; v += 2) {

	    if ((i == 0) && (u == -1)) continue;
	    if ((j == 0) && (v == -1)) continue;
	    if ((i == j) && (u != v)) continue;

	    Pull2Set[m][0] = (n == 0) ? 0 : NumberOfMonomers - 1;
	    Pull2Set[m][1] = (n == 0) ? 1 : -1;
	    Pull2Set[m][2] = i;
	    Pull2Set[m][3] = u;
	    Pull2Set[m][4] = j;
	    Pull2Set[m][5] = v;

// 	    printf("set %3ld : %3d %3d %3d %3d %3d %3d\n",
// 		   m,
// 		   Pull2Set[m][0],
// 		   Pull2Set[m][1],
// 		   Pull2Set[m][2],
// 		   Pull2Set[m][3],
// 		   Pull2Set[m][4],
// 		   Pull2Set[m][5]);

	    m++;

	  }
	}
      }
    }
  }

  if (m != MaxPull2Moves) ;

  printf("done\n");
  printf("--> max. number of pull-2 moves = %ld\n\n", MaxPull2Moves);

}

  // does various model integrity checks (in particular, checks the
  // consistency between the two representations of storing the
  // protein/polymer configuration, i.e. d-dimensional vectors for the
  // monomer coordinates vs entries in a d-dimensional occupancy
  // field)
  void CheckModelIntegrity()
{

  int m, n;

  for (n = 0; n < NumberOfMonomers - 1; n++)
    for (m = n + 1; m < NumberOfMonomers; m++)
      if (coords[n] -> Overlap(*coords[m])) ;

  m = 0;
  for (n = 0; n < NumberOfMonomers - 1; n++) {
    if (coords[n] -> Distance2(*coords[n + 1]) != 1);
    if ((seq[n] == hydrophobic) &&
	(seq[n + 1] == hydrophobic)) m++;
  }

  if (InternalHHContacts != m) ;

  m = 0;
  for (n = 0; n < NumberOfMonomers; n++) {
    if (site[n] != GetSite(*coords[n])) ;
    if (OccupancyField -> CheckKey(site[n]) != n) ;
    if ((seq[n] != hydrophobic) && (seq[n] != polar)) ;
    if (seq[n] == hydrophobic) m++;
  }

  if (NumberOfHMonomers != m) ;

  if (OccupancyField -> GetFillSize() != (unsigned long int)(NumberOfMonomers)) ;

  printf("\nModel integrity ok!\n");

}
  // returns the monomer location in the d-dimensional occupancy field
  // from the d-dimensional vector for the monomer coordinates
  // 1. reference to a "Vector" object (monomer coordinates)
  unsigned long int GetSite(Vector& v)
{

  int i;
  long int s;
  unsigned long int index = 0;

  for (i = 0; i < PolymerDim; i++) {
    s = v.Elem(i) % LSize;
    index += LOffset[i] * (unsigned long int)((s >= 0) ? s : s + LSize);
  }

  return index;

}


  // performs a pull move (labeled a "pull-2" move because there are
  // generally two monomers which are pulled at a time)
  void Pull2Move()
{

  int n, d, i, u, j, v;
  int dim, c, k, end, flag = 1;
  CoordsType dir;
  Vector* t;
  unsigned long int s;

  long int m = (long int)(gsl_rng_uniform(rng) * MaxPull2Moves);

  n = Pull2Set[m][0];
  d = Pull2Set[m][1];
  i = Pull2Set[m][2];
  u = Pull2Set[m][3];
  j = Pull2Set[m][4];

  // internal pull-2 move
  if (j == -1) {

    dim = coords[n] -> GetNormalDim(*coords[n - d]);
    i = (dim + i) % PolymerDim;

    // site L
    tmpcoords[n] -> Copy(*coords[n - d]);
    tmpcoords[n] -> Elem(i) += u;
    tmpsite[n] = GetSite(*tmpcoords[n]);
    if (OccupancyField -> Query(tmpsite[n]) != -1) return;

    // site C
    k = n + d;
    tmpcoords[k] -> Copy(*coords[n]);
    tmpcoords[k] -> Elem(i) += u;
    tmpsite[k] = GetSite(*tmpcoords[k]);
    c = OccupancyField -> Query(tmpsite[k]);
    if (c != -1) {
      if (c != k) return;
      else flag = 0;
    }

  }

  // end pull-2 move
  else {

    v = Pull2Set[m][5];
    k = n + d;
    dim = coords[n] -> GetNormalDim(*coords[k]);
    dir = coords[n] -> Elem(dim) - coords[k] -> Elem(dim);

    // site C
    tmpcoords[k] -> Copy(*coords[n]);
    tmpcoords[k] -> Elem((dim + i) % PolymerDim) += u * dir;
    tmpsite[k] = GetSite(*tmpcoords[k]);
    if (OccupancyField -> Query(tmpsite[k]) != -1) return;

    // site L
    tmpcoords[n] -> Copy(*tmpcoords[k]);
    tmpcoords[n] -> Elem((dim + j) % PolymerDim) += v * dir;
    tmpsite[n] = GetSite(*tmpcoords[n]);
    if (OccupancyField -> Query(tmpsite[n]) != -1) return;

  }

  BeadsIncr = d;
  BeadsStart = n;

  if (flag) {

    if (d == 1) end = NumberOfMonomers - 1;
    else end = 0;

    n += d;
    while ((n != end) && (tmpcoords[n] -> Distance2(*coords[n + d]) != 1)) {
      n += d;
      tmpsite[n] = site[n - 2*d];
      tmpcoords[n] -> Copy(*coords[n - 2*d]);
    }

  }

  BeadsEnd = n;
  BeadsInterval = BeadsIncr * (BeadsEnd - BeadsStart) + 1;

  // change of observables from "old" polymer fragment...
  if (BeadsInterval < BeadsIntervalThreshold) {
    if (H0PModelFlag)
      GetObservablesPartH0P(BeadsStart, BeadsEnd, -1);
  }

  // apply pull-2 move...
  n = BeadsStart;

  while (1) {

    OccupancyField -> Delete(site[n]);
    s = site[n];
    site[n] = tmpsite[n];
    tmpsite[n] = s;
    OccupancyField -> Insert(site[n], n);

    t = coords[n];
    coords[n] = tmpcoords[n];
    tmpcoords[n] = t;

    if (n == BeadsEnd) break;

    n += BeadsIncr;

  }

  // change of observables from "new" polymer fragment...
  if (BeadsInterval < BeadsIntervalThreshold) {
    if (H0PModelFlag)
      GetObservablesPartH0P(BeadsStart, BeadsEnd, +1);
  }
  else {
    if (H0PModelFlag)
      GetObservablesH0P();
  }

  MoveProposal = 1;

}
  // performs a bond-rebridging move
  void BondRebridgingMove()
{

  int i, n;
  int c1(-1), c2(-1), c3(-1), c4(-1), cdir(-1);
  int j1(-1), j2(-1), j3(-1), j4(-1), jdir(-1);
  int cmin, cmax;
  int ndim, d, dim, dir;
  Vector** tt;
  unsigned long int* s;

  n = (int)(gsl_rng_uniform(rng) * (NumberOfMonomers + 1));

  if (n > 1) {  // chain internal bond-rebridging

    c1 = n - 2;
    c2 = n - 1;

    ndim = coords[c1] -> GetNormalDim(*coords[c2]);
    d = (int)(gsl_rng_uniform(rng) * (2 * PolymerDim - 2));
    dim = d / 2;
    if (dim >= ndim) dim++;
    dir = 1 - 2 * (d % 2);

    coords[c1] -> Elem(dim) += dir;
    c3 = OccupancyField -> Query(GetSite(*coords[c1]));
    coords[c1] -> Elem(dim) -= dir;
    if (c3 == -1) return;

    coords[c2] -> Elem(dim) += dir;
    c4 = OccupancyField -> Query(GetSite(*coords[c2]));
    coords[c2] -> Elem(dim) -= dir;
    if (c4 == -1) return;

    if (abs(c4 - c3) != 1) return;

    if ((c4 - c3) == -1) {  // anti-parallel strands

      cdir = 1;

      if (c1 < c3) {
	cmin = c2;
	cmax = c4;
      }
      else {
	cmin = c3;
	cmax = c1;
      }

      if ((cmax - cmin) != 1) {
	j1 = cmin + (int)(gsl_rng_uniform(rng) * (cmax - cmin + 1));
	if (j1 != cmax) j2 = j1 + 1;
	else {
	  j2 = cmin;
	  if (cmin == c2) {
	    c2 = -1;
	    c4 = -1;
	  }
	  else {
	    c3 = -1;
	    c1 = -1;
	  }
	}
      }
      else {
	j1 = cmin;
	j2 = cmax;
      }

      ndim = coords[j1] -> GetNormalDim(*coords[j2]);
      d = (int)(gsl_rng_uniform(rng) * (2 * PolymerDim - 2));
      dim = d / 2;
      if (dim >= ndim) dim++;
      dir = 1 - 2 * (d % 2);

      coords[j1] -> Elem(dim) += dir;
      j3 = OccupancyField -> Query(GetSite(*coords[j1]));
      coords[j1] -> Elem(dim) -= dir;
      if ((j3 == -1) || ((cmin <= j3) && (j3 <= cmax))) return;

      coords[j2] -> Elem(dim) += dir;
      j4 = OccupancyField -> Query(GetSite(*coords[j2]));
      coords[j2] -> Elem(dim) -= dir;
      if ((j4 == -1) || ((cmin <= j4) && (j4 <= cmax))) return;

      if (abs(j4 - j3) != 1) return;

      if ((j4 - j3) == -1) jdir = 1;  // anti-parallel strands

    }

  }

  else {  // chain terminals bond-rebridging

    if (n == 0) {
      c1 = 0;
      c3 = 1;
    }
    else {
      c1 = NumberOfMonomers - 1;
      c3 = NumberOfMonomers - 2;
    }

    ndim = coords[c1] -> GetNormalDim(*coords[c3]);
    d = (int)(gsl_rng_uniform(rng) * (2 * PolymerDim - 1));
    if (d != (2 * PolymerDim - 2)) {
      dim = d / 2;
      if (dim >= ndim) dim++;
      dir = 1 - 2 * (d % 2);
    }
    else {
      dim = ndim;
      dir = coords[c1] -> Elem(ndim) - coords[c3] -> Elem(ndim);
    }

    coords[c1] -> Elem(dim) += dir;
    c3 = OccupancyField -> Query(GetSite(*coords[c1]));
    coords[c1] -> Elem(dim) -= dir;
    if (c3 == -1) return;

  }

  if (n == 0) {
    i = c3 - 1;
    d = -1;
  }
  else {
    i = 0;
    d = 1;
  }

  //  printf("c1-4, cdir : %3d %3d %3d %3d %3d\n", c1, c2, c3, c4, cdir);
  //  printf("j1-4, jdir : %3d %3d %3d %3d %3d\n", j1, j2, j3, j4, jdir);

  // apply bond-rebridging move...
  BeadsInterval = 0;
  for (n = 0; n < NumberOfMonomers; n++) {

    CoordsIndex[n] = i;
    CoordsPointer[n] = coords[i];
    tmpsite[n] = site[i];

    //    printf("Relabel: %3d --> %3d\n", n, i);

    if (n != i) {
      OccupancyField -> Replace(site[i], n);
      BeadsInterval++;
    }

    if (i == c1) {
      i = c3;
      c3 = -1;
      d *= cdir;
    }
    else if (i == c2) {
      i = c4;
      c4 = -1;
      d *= cdir;
    }
    else if (i == c3) {
      i = c1;
      c1 = -1;
      d *= cdir;
    }
    else if (i == c4) {
      i = c2;
      c2 = -1;
      d *= cdir;
    }
    else if (i == j1) {
      i = j3;
      j3 = -1;
      d *= jdir;
    }
    else if (i == j2) {
      i = j4;
      j4 = -1;
      d *= jdir;
    }
    else if (i == j3) {
      i = j1;
      j1 = -1;
      d *= jdir;
    }
    else if (i == j4) {
      i = j2;
      j2 = -1;
      d *= jdir;
    }
    else i += d;

  }

  tt = coords;
  coords = CoordsPointer;
  CoordsPointer = tt;

  s = site;
  site = tmpsite;
  tmpsite = s;

  // get "new" observables...
  // note: bond-rebridging moves do NOT alter
  //       the number of HH contacts of ISAWs
  if (H0PModelFlag) GetObservablesH0P();

  MoveProposal = 1;

}
  // performs a pivot move
  void PivotMove()
{

  int m, n, pp;
  long int po;
  Vector* t;
  unsigned long int s;

  // identity operation is excluded
  // because it does not generate a "new" configuration
  po = (long int)(gsl_rng_uniform(rng) * (NPiOps - 1));
  pp = (int)(gsl_rng_uniform(rng) * (NumberOfMonomers - 2)) + 1;

  if (pp < (NumberOfMonomers - 1 - pp)) {
    BeadsStart = pp - 1;
    BeadsIncr = -1;
    BeadsEnd = 0;
  }
  else {
    BeadsStart = pp + 1;
    BeadsIncr = +1;
    BeadsEnd = NumberOfMonomers - 1;
  }

  n = BeadsStart;
  while (1) {
    pivot[po] -> PivotOperation(*coords[pp], *coords[n], *tmpcoords[n]);
    tmpsite[n] = GetSite(*tmpcoords[n]);
    m = OccupancyField -> Query(tmpsite[n]);
    // if (m == "occupied") -->
    //  if (BeadsIncr == -1) AND (m > "pivot-point") --> return
    //    OR
    //  if (BeadsIncr == +1) AND (m < "pivot-point") --> return
    if (m != -1) {
      if (BeadsIncr == -1) {
	if (m > pp) return;
      }
      else if (m < pp) return;
    }
    if (n == BeadsEnd) break;
    n += BeadsIncr;
  }

  BeadsInterval = BeadsIncr * (BeadsEnd - BeadsStart) + 1;

  // change of observables from "old" polymer fragment...
  if (BeadsInterval < BeadsIntervalThreshold) {
    if (H0PModelFlag)
      GetObservablesPartH0P(BeadsStart, BeadsEnd, -1);
  }

  // apply pivot move...
  n = BeadsStart;
  while (1) {
    OccupancyField -> Delete(site[n]);
    t = coords[n];
    coords[n] = tmpcoords[n];
    tmpcoords[n] = t;
    if (n == BeadsEnd) break;
    n += BeadsIncr;
  }

  n = BeadsStart;
  while (1) {
    s = site[n];
    site[n] = tmpsite[n];
    tmpsite[n] = s;
    OccupancyField -> Insert(site[n], n);
    if (n == BeadsEnd) break;
    n += BeadsIncr;
  }

  // change of observables from "new" polymer fragment...
  if (BeadsInterval < BeadsIntervalThreshold) {
    if (H0PModelFlag)
      GetObservablesPartH0P(BeadsStart, BeadsEnd, +1);
  }
  else {
    if (H0PModelFlag)
      GetObservablesH0P();
  }

  MoveProposal = 1;

}

  // evaluates the number of non-bonded HH contacts (principal
  // observable) for the H0P model running over the entire protein
  // chain
  void GetObservablesH0P()
{

  int i, n, nn, nHH, nH0;
  long int HHContacts;
  long int H0Contacts;

  HHContacts = 0;
  H0Contacts = 0;

  for (n = 0; n < NumberOfMonomers; n++) {

    if (seq[n] == hydrophobic) {

      for (i = 0; i < PolymerDim; i++) {

	if ((coords[n] -> Elem(i) % LSize) != 0)
	  nn = OccupancyField -> Query(site[n] - LOffset[i]);
	else
	  nn = OccupancyField -> Query(site[n] + (LSize - 1) * LOffset[i]);

	if ((nn != -1) && (seq[nn] == hydrophobic)) HHContacts++;
  if ((nn != -1) && (seq[nn] == neutral)) H0Contacts++;
      }

    }

    if (seq[n] == neutral) {

      for (i = 0; i < PolymerDim; i++) {

	if ((coords[n] -> Elem(i) % LSize) != 0)
	  nn = OccupancyField -> Query(site[n] - LOffset[i]);
	else
	  nn = OccupancyField -> Query(site[n] + (LSize - 1) * LOffset[i]);

	if ((nn != -1) && (seq[nn] == hydrophobic)) H0Contacts++;

      }

    }

  }

  nHH = HHContacts - InternalHHContacts;
  nH0 = H0Contacts - InternalH0Contacts;
  observable[0] = 2 * nHH + nH0;

}

  // evaluates the number of non-bonded HH contacts (principal
  // observable) for the H0P model running only over the protein chain
  // fragment whose monomer locations have changed after a Monte Carlo
  // trial move
  void GetObservablesPartH0P(int nmin, int nmax, long int sign)
{

  int i, n, nn, nhH, nh0;
  long int HHint, HHext, H0int, H0ext;

  // in order to get the correct number of HH contacts for the
  // evaluated chain fragment, one needs to distinguish between
  // "internal" (i.e. within the chain fragment) and "external"
  // (i.e. with monomers outside of the chain fragment) HH contacts

  // --> "internal" HH contacts are counted twice (divide by 2)

  // note that "intra-chain" HH contacts cancel out when determining
  // the relative change of HH contacts (i.e. difference of HH
  // contacts between the "new" and "old" chain fragments)

  HHint = 0;  // "internal" HH contacts
  HHext = 0;  // "external" HH contacts
  H0int = 0;  // "internal" H0 contacts
  H0ext = 0;  // "eksternal" H0 contacts

  if (nmin > nmax) {
    n = nmax;
    nmax = nmin;
    nmin = n;
  }
  else n = nmin;

  for (; n <= nmax; n++) {

    if (seq[n] == hydrophobic) {

      for (i = 0; i < PolymerDim; i++) {

	      // check positive direction of coordinate...
	      if (((coords[n] -> Elem(i) + 1) % LSize) != 0)
	        nn = OccupancyField -> Query(site[n] + LOffset[i]);
	      else
	        nn = OccupancyField -> Query(site[n] - (LSize - 1) * LOffset[i]);

	      if ((nn != -1) && (seq[nn] == hydrophobic)) {
	        if ((nn < nmin) || (nn > nmax)) HHext++;
	        else HHint++;
        }
        if ((nn != -1) && (seq[nn] == neutral)) {
	        if ((nn < nmin) || (nn > nmax)) H0ext++;
	        else H0int++;
	      }

	      // check negative direction of coordinate...
	      if ((coords[n] -> Elem(i) % LSize) != 0)
	        nn = OccupancyField -> Query(site[n] - LOffset[i]);
	      else
	        nn = OccupancyField -> Query(site[n] + (LSize - 1) * LOffset[i]);

	      if ((nn != -1) && (seq[nn] == hydrophobic)) {
	        if ((nn < nmin) || (nn > nmax)) HHext++;
	        else HHint++;
        }
        if ((nn != -1) && (seq[nn] == neutral)) {
	        if ((nn < nmin) || (nn > nmax)) H0ext++;
	        else H0int++;
	      }

      }

    }

  }

  for (; n <= nmax; n++) {

    if (seq[n] == neutral) {

      for (i = 0; i < PolymerDim; i++) {

	      // check positive direction of coordinate...
	      if (((coords[n] -> Elem(i) + 1) % LSize) != 0)
	        nn = OccupancyField -> Query(site[n] + LOffset[i]);
	      else
	        nn = OccupancyField -> Query(site[n] - (LSize - 1) * LOffset[i]);

	      if ((nn != -1) && (seq[nn] == hydrophobic)) {
	        if ((nn < nmin) || (nn > nmax)) H0ext++;
	        else H0int++;
	      }

	      // check negative direction of coordinate...
	      if ((coords[n] -> Elem(i) % LSize) != 0)
	        nn = OccupancyField -> Query(site[n] - LOffset[i]);
	      else
	        nn = OccupancyField -> Query(site[n] + (LSize - 1) * LOffset[i]);

	      if ((nn != -1) && (seq[nn] == hydrophobic)) {
	        if ((nn < nmin) || (nn > nmax)) H0ext++;
	        else H0int++;
	      }

      }

    }

  }

  nhH += sign * (HHint / 2 + HHext);
  nh0 += sign * (H0int / 2 + H0ext);
  observable[0] = 2 * nhH + nh0;
}



};

class MonteCarlo {
/*
  MonteCarlo class:
  This class provides the important control parameters  and routines used for
  Wang-Landau sampling (e.g. initial/final modification factors, flatness criterion etc.) 
*/

public:
  // inisialisasi parameter-parameter dalam algoritma Wang-Landau
  int RngType = 1;    // tipe random number generator
  unsigned long int RngSeed = 10;      // benih random number generator                 

  double ModFactorInit = 1.0;                // initial modification factor (ln)
  double ModFactorDivider = 2.0;             // modification factor divider
  double ModFactorThreshold = 1.0e-8;        // final modification factor (ln)
  double FlatnessMeasure = 0.8;              // histogram flatness criterion
  unsigned long int MinHitsBin = 0;          // min. hits per histogram bin
  int HistogramCheckType = 0;                // histogram check type, see "Histogram" class
  unsigned long int HistogramCheckInterval = 1000000;   // # MC steps between successive flatness checks
  int DOSUpdateType = 1;      // DOS update type, see "Histogram" class
  unsigned long int DOSDumpInterval = 0;    // # MC steps between successive storage of "Histogram" data
  int MaskAllFlag = 0;        // flag: 1 = mask all histogram entries as visited

  unsigned long int MCMoves = 0;            // current # MC steps
  unsigned long int MCMovesAccepted = 0;    // current # accepted MC steps
  unsigned long int MCMovesMem = 0;     
  unsigned long int MCMovesIteration = 0;
  double ModFactor = ModFactorInit;         // current modification factor
  int WLIteration = 1;                      // current Wang-Landau iteration
  ObservableType Emin = 0;                  // current energy minimum


  // constructor
  MonteCarlo() {

  // inisialisasi random number generator
  switch (RngType) {
  case  1 : {rng = gsl_rng_alloc(gsl_rng_ranlxd1); break;}
  case  2 : {rng = gsl_rng_alloc(gsl_rng_mt19937); break;}
  default : {rng = gsl_rng_alloc(gsl_rng_ranlxd2);}
  }
  
  if (RngSeed == 0) RngSeed = 1;
  gsl_rng_set(rng, RngSeed);

}

  // destructor
  ~MonteCarlo(){
  gsl_rng_free(rng);    // "deletes" the GSL random number generator
}

  // Wang-Landau sampling (WLS)
  // 1. pointer to a "Histogram" object
  // 2. pointer to a "Model" object
  void WangLandauSampling(Histogram* h, H0PModel* m)
{

  unsigned long int n;
  long int prev, cur;

  m -> PrintObservables();
  prev = h -> GetIndex(m -> observable);

  while (ModFactor >= ModFactorThreshold) {

    for (n = 0; n < HistogramCheckInterval; n++) {

      MCMoves++;
      m -> DoMCMove();

      if (m -> MoveProposal) {

	cur = h -> GetIndex(m -> observable);

	if ((cur >= 0) && (gsl_rng_uniform(rng) < (exp(h -> GetDOS(prev) - h -> GetDOS(cur))))) {

	  prev = cur;
	  MCMovesAccepted++;

	  // in the "H0PModel" class the measured observable is "number
	  // of non-bonded HH contacts"; the energy is minus this
	  // quantity, therefore the minus sign
	  // --> needs to be adapted for each physical model
	  if (-m -> observable[0] < Emin) {
	    Emin = (-m -> observable[0]);
	    printf("New Emin = %ld ( MC moves = %lu )\n", Emin, MCMoves);
	    fflush(stdout);
	    m -> WriteState(1, "confsEmin.mol2");
	    m -> WriteState(3, "confsEmin.xyz");
	  }

	  h -> CheckItinerancy(prev, MCMoves, MCMovesMem);

	}

	else m -> UnDoMCMove();
  
  
  if (WLIteration >= 26){
    m -> GetRadiusofGyration();
    m -> GetTortuosity();
  }
  

      }
      
      h -> Update(prev, n, DOSUpdateType, ModFactor);

      if ((DOSDumpInterval) && ((MCMoves % DOSDumpInterval) == 0))
	h -> SaveState(MCMoves / DOSDumpInterval, NULL, false);

    }

    if (h -> CheckHistogram(HistogramCheckType, FlatnessMeasure, MinHitsBin)) {
      printf("DOS: iteration = %d ( mod. factor = %15.8e , MC moves : %lu %lu )\n",
	     WLIteration, ModFactor, MCMoves, MCMoves - MCMovesIteration);
      h -> SaveState(WLIteration, "hdata_iteration");
      h -> PrintNormDOS();
      WLIteration++;
      MCMovesIteration = MCMoves;
      ModFactor /= ModFactorDivider;
      h -> ResetHist();
    }
    //    else printf("Histogram not flat ( MC moves = %lu )\n", MCMoves);

    h -> SaveState();
    m -> WriteState(4, "coords_current.xyz");
    fflush(stdout);

  }

  printf("Wang-Landau (WL) sampling done.\n");
  printf("MC moves = %lu\n", MCMoves);
  printf("Accepted MC moves = %lu\n", MCMovesAccepted);
  printf("Acceptance ratio = %15.8e\n", double(MCMovesAccepted) / double(MCMoves));
  
  
  h -> SaveState(0, "hdata_final.dat");
  h -> PrintNormDOS("dos.dat");
  
  
  if (WLIteration >= 26){
    m -> GetAverageRg();
    m -> GetAverageTau();
  }
  m -> printRg();
  m -> printTau();

}

};



// this is the main program
int main(int argc, char* argv[])
{

  MonteCarlo* mcs;
  Histogram* histogram = NULL;
  H0PModel* model;

  // instantiate a "MonteCarlo" class object
  mcs = new MonteCarlo();

  // instantiate a "physical model" class object
  model = new H0PModel(argv[1]);

  // instantiate a "Histogram" class object for Wang-Landau sampling
  histogram = new Histogram();
  
  // run Monte Carlo simulation, either Wang-Landau sampling (default)
  mcs -> WangLandauSampling(histogram, model);

  delete model;
  if (histogram != NULL) delete histogram;
  delete mcs;

  printf("\nSimulation successfully completed.\n");

}
