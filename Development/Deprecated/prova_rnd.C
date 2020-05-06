/*
 * =====================================================================================
 *
 *       Filename:  prova_rnd.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/04/2020 17:12:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Luigi Berducci (), berducci.1647887@studenti.uniroma1.it
 *        Company:  University of Rome La Sapienza
 *
 * =====================================================================================
 */
#include <TStopwatch.h>

void prova_rnd(){
    TStopwatch timer;
    TRandom rnd;
    const long int n=1000000;
    array<double, n> A;
    for(long int i=0; i<n; i++){
        A[i] = rnd.Gaus(.5, .5);
    }
    std::cout << "Time: " << timer.CpuTime() << std::endl;
}
