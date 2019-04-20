# Description 

This is the program to compute the Harborth constant g(G) by brute force, that accompanies our article: "On the Harborth constant of `C_3 ⊕ C_3n`

# Program requirements

Please note that the program has been tested using the following elements:

1- The used system is `ubuntu 16.04.10` with `kernel 4.15.0-39`
2- The used SHELL is `/bin/bash` version `4.3.48`
3- The `gcc` compiler is used with `5.4.0-6 version`

# Sample to run the program

1- Compile the program first `make -f makefile`
2- run the program with the compiled file `./sommenulle`

# The user manual

For now, There is no interface but you only have to change this source code to get results for the Harborth constant for the different finite abelian groups.

You have to change the parameters of: `MODULO1, MODULO2, MODULO3 and MODULO4`, according to the finite abelian group of which you are looking for its Harborth constant. 

## Example: 
For The group `C_3 ⊕ C_12`, you have to enter the following parameters:
 
## Initialization: 
```
define RANG 4

define MODULO1 1
define MODULO2 1
	if (RANG>=3)
		define MODULO3 3
	endif
	if (RANG==4)
		define MODULO4 12
	endif
	if (RANG==2)
		define CARDINAL (MODULO1*MODULO2)
		define EXPOSANT MODULO2
	endif
	if (RANG==3)
		define CARDINAL (MODULO1*MODULO2*MODULO3)
		define EXPOSANT MODULO3
	endif
	if(RANG==4)
		define CARDINAL (MODULO1*MODULO2*MODULO3*MODULO4)
		define EXPOSANT MODULO4
	endif
```
