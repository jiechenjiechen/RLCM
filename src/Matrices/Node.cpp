#include "Node.hpp"

#define INITVAL_NumChild  0
#define INITVAL_n         0
#define INITVAL_r         0
#define INITVAL_start     -1
#define INITVAL_LogAbsDet 0.0
#define INITVAL_Sign      1
#define INITVAL_offset    DBL_MAX
#define INITVAL_Dim       0
#define INITVAL_PartDim   -1
#define ADDRESS_WIDTH     14
#define INTERVAL_WIDTH    7


//--------------------------------------------------------------------------
Node::
Node() {
  Parent       = NULL;
  LeftChild    = NULL;
  RightSibling = NULL;
  pivots       = NULL;
  BBox         = NULL;
  NumChild     = INITVAL_NumChild;
  n            = INITVAL_n;
  r            = INITVAL_r;
  start        = INITVAL_start;
  offset       = INITVAL_offset;
  Dim          = INITVAL_Dim;
  PartDim      = INITVAL_PartDim;
  mLogDet      = { INITVAL_LogAbsDet, INITVAL_Sign };
  Init();
}


//--------------------------------------------------------------------------
void Node::
Init(void) {
  ReleaseAllMemory();
  NumChild  = INITVAL_NumChild;
  n         = INITVAL_n;
  r         = INITVAL_r;
  start     = INITVAL_start;
  offset    = INITVAL_offset;
  Dim       = INITVAL_Dim;
  PartDim   = INITVAL_PartDim;
  mLogDet   = { INITVAL_LogAbsDet, INITVAL_Sign };
}


//--------------------------------------------------------------------------
void Node::
ReleaseAllMemory(void) {
  Delete_1D_Array<INTEGER>(&pivots);
  Delete_1D_Array<double>(&BBox);
  Dim = INITVAL_Dim;
  A.ReleaseAllMemory();
  Sigma.ReleaseAllMemory();
  W.ReleaseAllMemory();
  U.ReleaseAllMemory();
  Z.ReleaseAllMemory();
  V.ReleaseAllMemory();
  c.ReleaseAllMemory();
  d.ReleaseAllMemory();
  C.ReleaseAllMemory();
  D.ReleaseAllMemory();
  E.ReleaseAllMemory();
  Xi.ReleaseAllMemory();
  Theta.ReleaseAllMemory();
  normal.ReleaseAllMemory();
  FACT.ReleaseAllMemory();
  P.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
Node::
Node(const Node &G) {
  Parent       = NULL;
  LeftChild    = NULL;
  RightSibling = NULL;
  pivots       = NULL;
  BBox         = NULL;
  NumChild     = INITVAL_NumChild;
  n            = INITVAL_n;
  r            = INITVAL_r;
  start        = INITVAL_start;
  offset       = INITVAL_offset;
  Dim          = INITVAL_Dim;
  PartDim      = INITVAL_PartDim;
  mLogDet      = { INITVAL_LogAbsDet, INITVAL_Sign };
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
Node& Node::
operator= (const Node &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void Node::
DeepCopy(const Node &G) {

  if (NumChild != G.NumChild || r != G.r) {
    ReleaseAllMemory();
  }

  if (G.pivots) {
    if (!pivots) {
      New_1D_Array<INTEGER, INTEGER>(&pivots, G.r);
    }
    memcpy(pivots, G.pivots, G.r*sizeof(INTEGER));
  }
  else {
    Delete_1D_Array<INTEGER>(&pivots);
  }

  if (G.BBox) {
    if (!BBox) {
      New_1D_Array<double, INTEGER>(&BBox, G.Dim*2);
    }
    memcpy(BBox, G.BBox, G.Dim*2*sizeof(double));
  }
  else {
    Delete_1D_Array<double>(&BBox);
  }

  Parent       = G.Parent;
  NumChild     = G.NumChild;
  LeftChild    = G.LeftChild;
  RightSibling = G.RightSibling;
  n            = G.n;
  r            = G.r;
  start        = G.start;
  A            = G.A;
  Sigma        = G.Sigma;
  W            = G.W;
  U            = G.U;
  Z            = G.Z;
  V            = G.V;
  c            = G.c;
  d            = G.d;
  C            = G.C;
  D            = G.D;
  E            = G.E;
  Xi           = G.Xi;
  Theta        = G.Theta;
  normal       = G.normal;
  offset       = G.offset;
  FACT         = G.FACT;
  Dim          = G.Dim;
  PartDim      = G.PartDim;
  P            = G.P;
  mLogDet      = { G.mLogDet.LogAbsDet, G.mLogDet.Sign };

}


//--------------------------------------------------------------------------
Node::
~Node() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
void Node::
PrintNode(INTEGER MaxNumChild, INTEGER MaxN) const {

  std::string format;
  int len;
  char slen[1000];

  //---------- Address -----------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("Address").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("p  ");
  printf(format.c_str(), this);

  //---------- Parent ------------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("Parent").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("p  ");
  printf(format.c_str(), Parent);

  //---------- NumChild ----------------------------------------------------
  sprintf(slen, "%ld", (long)MaxNumChild);
  len = std::max<int>(std::string(slen).length(),
                      std::string("#Child").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("d  ");
  printf(format.c_str(), NumChild);

  //---------- LeftChild ---------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("LeftChild").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("p  ");
  printf(format.c_str(), LeftChild);

  //---------- RightSibling ------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("RightSibling").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("p  ");
  printf(format.c_str(), RightSibling);

  //---------- n -----------------------------------------------------------
  sprintf(slen, "%ld", (long)MaxN);
  len = std::max<int>(std::string(slen).length(),
                      std::string("n").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("d  ");
  printf(format.c_str(), n);

  //---------- start -------------------------------------------------------
  sprintf(slen, "%ld", (long)MaxN);
  len = std::max<int>(std::string(slen).length(),
                      std::string("start").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("d  ");
  printf(format.c_str(), start);

  //---------- end ---------------------------------------------------------
  sprintf(slen, "%ld", (long)MaxN);
  len = std::max<int>(std::string(slen).length(),
                      std::string("end").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("d");
  printf(format.c_str(), start+n-1);

  //---------- BBox --------------------------------------------------------
  if (BBox != NULL) {
    printf("  ");
    for (INTEGER i = 0; i < Dim; i++) {
      printf("[%+5.3f, %+5.3f]", BBox[i], BBox[i+Dim]);
      if (i != Dim-1) {
        printf(" x ");
      }
    }
  }

  //------------------------------------------------------------------------
  printf("\n");

}


//--------------------------------------------------------------------------
INTEGER Node::
GetMemEstNode(void) const {
  INTEGER mem = 0;
  if (NumChild == 0) {
    mem += Square(n); // A
    mem += n*r; // U
  }
  if (NumChild != 0) {
    mem += Square(r); // Sigma
  }
  if (NumChild != 0 && Parent != NULL) {
    mem += Square(r); // W
  }
  return mem;
}


//--------------------------------------------------------------------------
INTEGER Node::
GetMemEstTree(void) const {
  INTEGER i = 0;
  Node *child = LeftChild;
  INTEGER mem = 0;
  for (i = 0; i < NumChild; i++, child = child->RightSibling) {
    mem += child->GetMemEstTree();
  }
  mem += GetMemEstNode();
  return mem;
}


//--------------------------------------------------------------------------
INTEGER Node::
GetMemEstNodeAonly(void) const {
  INTEGER mem = 0;
  if (NumChild == 0) {
    mem = Square(n); // A
  }
  return mem;
}


//--------------------------------------------------------------------------
INTEGER Node::
GetMemEstTreeAonly(void) const {
  INTEGER i = 0;
  Node *child = LeftChild;
  INTEGER mem = 0;
  for (i = 0; i < NumChild; i++, child = child->RightSibling) {
    mem += child->GetMemEstTreeAonly();
  }
  mem += GetMemEstNodeAonly();
  return mem;
}


//--------------------------------------------------------------------------
void Node::
CopyTree(Node **mNode) {

  *mNode = new Node;
  (*mNode)->DeepCopy(*this);

  Node *child = LeftChild;
  Node *mNodeChild = NULL;
  Node *LastChild = NULL;
  if (NumChild != 0) {
    child->CopyTree(&mNodeChild);
    mNodeChild->Parent = *mNode;
    (*mNode)->LeftChild = mNodeChild;
    LastChild = mNodeChild;
    child = child->RightSibling;
  }

  for (INTEGER i = 1; i < NumChild; i++, child = child->RightSibling) {
    child->CopyTree(&mNodeChild);
    mNodeChild->Parent = (*mNode);
    LastChild->RightSibling = mNodeChild;
    LastChild = mNodeChild;
  }

}


//--------------------------------------------------------------------------
void Node::
DestroyTree(void) {
  if (RightSibling) {
    RightSibling->DestroyTree();
    delete RightSibling;
    RightSibling = NULL;
  }
  if (LeftChild) {
    LeftChild->DestroyTree();
    delete LeftChild;
    LeftChild = NULL;
  }
}


//--------------------------------------------------------------------------
void Node::
TreeReleaseAllMemory(void) {
  if (RightSibling) {
    RightSibling->TreeReleaseAllMemory();
  }
  if (LeftChild) {
    LeftChild->TreeReleaseAllMemory();
  }
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
void Node::
PrintTree(INTEGER MaxNumChild, INTEGER MaxN) const {

  // Print the field names
  std::string format;
  int len;
  char slen[1000];

  //---------- Address -----------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("Address").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s  ");
  printf(format.c_str(), std::string("Address").c_str());

  //---------- Parent ------------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("Parent").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s  ");
  printf(format.c_str(), std::string("Parent").c_str());

  //---------- NumChild ----------------------------------------------------
  sprintf(slen, "%ld", (long)MaxNumChild);
  len = std::max<int>(std::string(slen).length(),
                      std::string("#Child").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s  ");
  printf(format.c_str(), std::string("#Child").c_str());

  //---------- LeftChild ---------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("LeftChild").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s  ");
  printf(format.c_str(), std::string("LeftChild").c_str());

  //---------- RightSibling ------------------------------------------------
  len = std::max<int>(ADDRESS_WIDTH, std::string("RightSibling").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s  ");
  printf(format.c_str(), std::string("RightSibling").c_str());

  //---------- n -----------------------------------------------------------
  sprintf(slen, "%ld", (long)MaxN);
  len = std::max<int>(std::string(slen).length(),
                      std::string("n").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s  ");
  printf(format.c_str(), std::string("n").c_str());

  //---------- start -------------------------------------------------------
  sprintf(slen, "%ld", (long)MaxN);
  len = std::max<int>(std::string(slen).length(),
                      std::string("start").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s  ");
  printf(format.c_str(), std::string("start").c_str());

  //---------- end ---------------------------------------------------------
  sprintf(slen, "%ld", (long)MaxN);
  len = std::max<int>(std::string(slen).length(),
                      std::string("end").length());
  sprintf(slen, "%d", len);
  format = std::string("%") + std::string(slen) + std::string("s");
  printf(format.c_str(), std::string("end").c_str());

  //---------- BBox --------------------------------------------------------
  if (BBox != NULL) {
    printf("  BBox");
  }

  //------------------------------------------------------------------------
  printf("\n");

  // Recursively print node information
  PrintTreeDownward(MaxNumChild, MaxN);

}


//--------------------------------------------------------------------------
void Node::
PrintTreeDownward(INTEGER MaxNumChild, INTEGER MaxN) const {
  PrintNode(MaxNumChild, MaxN);
  INTEGER i = 0;
  Node *child = LeftChild;
  for (i = 0; i < NumChild; i++, child = child->RightSibling) {
    child->PrintTreeDownward(MaxNumChild, MaxN);
  }
}
