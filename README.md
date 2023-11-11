# PLC Project 1
Write a program in the Julia programming language to locate an object in three-dimensional space using the time difference of a signal collected between an object and four satellites.

# Input / Output

All the input is US-ASCII text (decimal fractions, actually) to be read from the standard input stream. In the input are decimal factions separated by white space. All the output is to be written to the standard output stream.

The input consists of two parts: the location of four satellites, followed by several lines of times. The locations of the satellites in three-dimensional Euclidean space are in meters.

```
 -29891879.65      29967333.82     -25251373.37
 -28496603.73     -30010958.93     -26741603.64
 -22583800.11     -28474887.66      30378930.53
 -26176414.09      33705464.57      22747328.98
 194930775.065889 174716098.262710 150139783.885383 180201150.009913
 191784661.765340 172603245.743806 146792018.829073 176002731.614972
```

The location of the four satellites remains fixed for all of the input. Each additional line after the first four represents some object that you are to locate. These lines containing four times. These times are not the time it takes the signal to reach the object. These times are timestamps in nanosecodnds noted on arrival in some common frame of reference. Adding a constant to each of the four timestamps does not change the answer.

The output consists of a few of the intermediate values (g, h, j, m, and o). These values are a check on the calculation and should be printed in scientific notation with just two digits to the right of the decimal. Because of the quadratic formula there are two possible answers. Print both locations. Print x, y, and z to the nearest meter. Also print r -- the distance to the origin.

```
g= -9.52e+02, h=  4.01e+09, j= -2.29e+08, m=  2.55e+20, o=  4.46e+33
+) x=    5968787, y=   -5552597, z=    4200921; r=    9170905
-) x=   46595294, y=   -7819528, z=    4158246; r=   47429500

g=  2.98e+02, h= -1.35e+09, j=  6.46e+07, m=  2.55e+19, o=  5.47e+32
+) x=   47520917, y=   -7453769, z=    4701171; r=   48331121
-) x=    4430167, y=   -5235480, z=    4556385; r=    8233910
```

The data is so constructed such that the correct position can be caluculated using the IEEE 754 64-bit floating-point data type. 