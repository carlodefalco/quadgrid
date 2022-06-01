# quadgrid
simple cartesian quad grid for c++/octave

`quadgrid_cpp.h`      contiene la classe (template) `quadgrid_t`
`quadgrid_cpp.cpp`    è una demo di come utilizzare la classe in una applicazione c++

`quadgrid.h`          definisce la classe `quadgrid` derivata da `octave_base_value` che permette di definire un oggetto `quadgrid_t` in Octave
`quadgrid.cc`         definisce le due funzioni Octave `quadgrid` (costruttore della classe quadgrid) e `quadgrid_loop` (esempio di uso di quadgrid in Octave)

per compilare le demo in c++

    mpicxx -std=c++17 -I. -o quadgrid_cpp quadgrid_cpp.cpp 
    mpicxx -std=c++17 -I. -o mass_matrix_example mass_matrix_example.cpp 
    
per compilare le funzioni utilizzabili da Octave

    CXX=mpicxx CPPFLAGS="-I. -std=c++17" mkoctfile quadgrid.cc 
    
per testare da Octave

     q = quadgrid (10, .1, 5, .02)
     quadgrid_loop (q)
     
