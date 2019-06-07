# BdG_cpp
This solves the Bogoliubov-de Gennes equations and gap equations in the s-wave superconductor with the use of the Reduced-Shifted Conjugate-Gradient Method method. See, Y. Nagai, Y. Shinohara, Y. Futamura, and T. Sakurai,[arXiv:1607.03992v2 or DOI:10.7566/JPSJ.86.014708]. http://dx.doi.org/10.7566/JPSJ.86.014708  

This code is for understanding the methods. 
The code was written by c++11 with Boost 1.69 and MKL.

So the compile command is like 

```
g++ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -lm -lpthread -liomp5  -I/usr/local/Cellar/boost/1.69.0_2/include testbdg.cpp
```

I do not guarantee any results calculated by this code. 
