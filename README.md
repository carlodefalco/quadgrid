# quadgrid
## simple cartesian quad grid with embedded particles for c++/octave

This library provides simple (template) classes that for a useful
starting point for the implementation of <a url="https://doi.org/10.1016/bs.aams.2020.12.002">MPM/PIC</a> methods 

The main 
* `quadgrid_cpp.h`            contiene la classe (template) `quadgrid_t`
* `particles.h`               contiene la classe `particles_t` che rappresenta un insieme di particelle distribuite nella griglia
* `quadgrid_cpp.cpp`          è una demo di come utilizzare la classe in una applicazione c++
* `mass_matrix_example.cpp`   è una demo di come utilizzare la classe in una applicazione c++ per costruire la matrice di massa
* `particle_sort_example.cpp` è una demo di come utilizzare la classe in una applicazione c++ per gestire particelle e interpolazione

* `quadgrid.h`          definisce la classe `quadgrid` derivata da `octave_base_value` che permette di definire un oggetto `quadgrid_t` in Octave
* `quadgrid.cc`         definisce le due funzioni Octave `quadgrid` (costruttore della classe quadgrid) e `quadgrid_loop` (esempio di uso di quadgrid in Octave)

per compilare le demo in c++

    mpicxx -std=c++17 -I. -o quadgrid_cpp quadgrid_cpp.cpp 
    mpicxx -std=c++17 -I. -o mass_matrix_example mass_matrix_example.cpp
    mpicxx -std=c++17 -I. -o particle_sort_example particle_sort_example.cpp 
    
per compilare le funzioni utilizzabili da Octave

    CXX=mpicxx CPPFLAGS="-I. -std=c++17" mkoctfile quadgrid.cc 
    
per testare da Octave

     q = quadgrid (10, .1, 5, .02)
     quadgrid_loop (q)
     
