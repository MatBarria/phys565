#include "../lib/CreateTuples.h"
#include <TStopwatch.h>

int main(int argc, char *argv[]) {

    if (argc != 5) {
        std::cerr << "Please give 4 arguments: input file, output directory, "
                     "channel, JES"
                  << std::endl;
        return -1;
    }

    TString input(argv[1]);
    TString output(argv[2]);
    TString channel(argv[3]);
    float JES = std::stof(argv[4]);

    std::cout << "Input: " << input << std::endl;
    std::cout << "Output: " << output << std::endl;
    std::cout << "Channel: " << channel << std::endl;
    std::cout << "JES: " << JES << std::endl;
    // std::cout << std::boolalpha;
    // std::cout << "Is data? " << is_data << std::endl;
    // std::cout << "Is signal? " << is_signal << std::endl;

    // TStopwatch timer = TStopwatch();
    // timer.Start();
    TStopwatch timer;
    timer.Start();

    CreateTuple create_tuple = CreateTuple(input, output, channel, JES);

    create_tuple.setBranchesAddressesInput();
    create_tuple.setBranchesAddressesOutput();
    create_tuple.fillOutputTree(channel);
    create_tuple.saveTree();
    timer.Print();
    return 0;
}
