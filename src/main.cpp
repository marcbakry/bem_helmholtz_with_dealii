#include <cstdlib>
#include <iostream>

#include "acoustic_bem.hpp"

int main()
{
    std::cout << "+------------------------------------------------+" << std::endl;
    std::cout << "| Resolution du probleme de Dirichlet acoustique |" << std::endl;
    std::cout << "| par methode BEM (resolution du simple couche)  |" << std::endl;
    std::cout << "+------------------------------------------------+" << std::endl;
    AcousticBEM myBEM(10.0);
    myBEM.run();
    std::cout << "+------------------------------------------------+" << std::endl;
    return EXIT_SUCCESS;
}