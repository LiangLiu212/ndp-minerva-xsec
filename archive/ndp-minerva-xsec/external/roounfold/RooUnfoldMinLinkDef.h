// ROOT dictionary request for the RooUnfold kBayes build. Every compiled class
// that carries a ClassDef needs a dictionary entry, or its Streamer() is an
// undefined symbol and the library fails to load (ROOT resolves at load time).
// We only ever instantiate RooUnfoldResponse + RooUnfoldBayes, but the base
// RooUnfold::New() pulls Svd/Invert/BinByBin into the vtable, so dict them too.
#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class RooUnfold+;
#pragma link C++ class RooUnfoldResponse+;
#pragma link C++ class RooUnfoldBayes+;
#pragma link C++ class RooUnfoldErrors+;
#pragma link C++ class RooUnfoldInvert+;
#pragma link C++ class RooUnfoldBinByBin+;
#pragma link C++ class RooUnfoldSvd+;
#endif
