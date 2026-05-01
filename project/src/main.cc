#include "../lib/CreateTuples.h"
#include <TStopwatch.h>

int main(int argc, char *argv[]) {

    if (argc != 4) {
        std::cerr << "Please give 3 arguments: input file, output file "
                  << std::endl;
        return -1;
    }

    TString input(argv[1]);
    TString output(argv[2]);
    //TString era(argv[3]);
    TString channel(argv[3]);
    //const bool is_data = *argv[5] == 'T';
    //const bool is_signal = *argv[6] == 'T';

    std::cout << "Input: " << input << std::endl;
    std::cout << "Output: " << output << std::endl;
    std::cout << "Channel: " << channel << std::endl;
    //std::cout << std::boolalpha;
    //std::cout << "Is data? " << is_data << std::endl;
    //std::cout << "Is signal? " << is_signal << std::endl;

    // TStopwatch timer = TStopwatch();
    // timer.Start();
    TStopwatch timer;
    timer.Start();

    CreateTuple create_tuple =
        CreateTuple(input, output, channel);

    create_tuple.setBranchesAddressesInput();
    create_tuple.setBranchesAddressesOutput();
    create_tuple.fillOutputTree(channel);
    create_tuple.saveTree();
    timer.Print();
    return 0;
}
