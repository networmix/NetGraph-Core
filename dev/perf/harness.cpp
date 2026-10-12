// The A/B/A timing harness: the rows of aba_rows.hpp compiled against one Core build's include/
// and src/ (build.sh). harness_main (aba_rows.hpp) reads the command line:
//
//   harness aba-list | aba-run <row> <iterations> | aba-calibrate <row> <seconds>
#include "aba_rows.hpp"

int main(int argc, char** argv) { return abarows::harness_main(argc, argv, abarows::rows()); }
