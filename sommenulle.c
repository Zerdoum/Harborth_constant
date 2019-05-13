// ************************* Introduction *************************

// This is the program to compute the Harborth constant g(G) by brute force, that accompanies our article: "On the Harborth constant of C_3 ⊕ C_3n".

// G is a finite abelian group.

// The Harborth constant g(G) is the smallest integer k such that each set over G whith size at least k, has a subset of size e=exp(G) that sums to 0. The Harborth constant g(G) corresponds to e= exp(G).

// There is no interface but you only have to change this source code to get results for the Harborth constant for the different finite abelian groups.

// This program is valid for any finite abelian group. With the hardware at our disposal it is possible to compute the Harborth constant for finite abelian groups of order up to about 45.

// The main limiting factor is memory.

// In order to increase the size of accessible groups,currently, we are working on a another version, more efficient based on data compression.

// The group intervenes only in this step [Initialization].

// The subsets of G are represented by a bitmap. 

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <malloc.h>
#include <time.h>

typedef long long unsigned u64;
typedef unsigned char u8;

// ************************* The user manual ************************* 

// [Initialization]: in this step you have to change the parameters of: MODULO1, MODULO2, MODULO3 and MODULO4, according to the finite abelian group of which you are looking for its Harborth constant. 

// Examples: For The group C_3 ⊕ C_12, you have to enter the following parameters:
 
//[Initialization].
#define RANG 2

#define MODULO1 3
#define MODULO2 12
#if (RANG>=3)
#define MODULO3 0
#endif
#if (RANG==4)
#define MODULO4 0
#endif
#if (RANG==2)
#define CARDINAL (MODULO1*MODULO2)
#define EXPOSANT MODULO2
#endif
#if (RANG==3)
#define CARDINAL (MODULO1*MODULO2*MODULO3)
#define EXPOSANT MODULO3
#endif
#if(RANG==4)
#define CARDINAL (MODULO1*MODULO2*MODULO3*MODULO4)
#define EXPOSANT MODULO4
#endif


// ************************* Functions ************************* 
// In this step you have no change to do. 

// Let's first look at the caclulation of w ---> x
// We consider a Pascal's triangle whose first column a(i,0)=1,
// The diagonal a(i,i)=2^i (powers of 2) and whose other values are defined by the Pascal's triangle relation. 
// a(i,j)=a(i-1,j)+a(i-1,j-1).
// Example :
// O
// 0 1
// 0 1 2
// 0 1 3 4
// 0 1 4 7 8 ...

// Let w be a bitmap, for example w=1100. The value x will be given by browsing this triangle.  
// The starting point is the position (-1,-1)
// We browse w from the left to the right.
// When we meet 1, we move diagonally down right. 
// When we meet 0, we move vertically down.
// The result is x which is the sum of all the summits encountered values.

// Example : w=1100
//-----------------
//
// O
//   \               1 : Move diagonally down right 
// 0  1
//      \            1 : Move diagonally down right
// 0  1  2
//       |           0 : Move vertically down
// 0  1  3  4
//       |           0 : Move vertically down 
// 0  1  4  7  8

// The result : the cumulation 0 + 1 + 2 + 3 + 4 = 10


// The function numsubset: return the number of the subset of a set of n elements. 
// The subsets of G are represented by a bitmap. 
u64 numsubset(u64 w,int n)
{
   	int i,j; // The index of Pascal's triangle. 
   	u64 a;   // Modified Pascal's triangle value. 
   	u64 p;   // Pascal's triangle value. 
   	u64 x;  // Cumulative of all the summits.
   	//[Initialization].
   	i=-1;
   	j=-1;
   	p=0;
   	a=0;
   	x=0;
   	while (n>0)
   	{
      	// According to the value of the number 
     		switch (w&1)
      		{
        	case 0:
          		++i;           // Move vertically down 
          		a=2*a-p;       // Update of the modified Pascal coefficient 
          		p=(p*i)/(i-j); // Update of Pascal coefficient
          		break;
        	case 1:
          		if (j==-1) p=1; else p=(p*(i-j))/(j+1);  // p(i-1,j) in terms of p(i-1,j-1)
          		a=2*a+p;        // Update of a
          		++i; ++j; // Diagonally down right
          		if (i==j) p=1; else p=(p*i)/(i-j);  // update of p(i,j) in terms of p(i-1,j)
          		break;
      		}
      		// shift w to position the next number  
      		w>>=1;
      		// cumulative value of the summit
      		x+=a;
      		--n; 
   	}
   	return x;
}


// The function subsetnum: return the subset of a given number 
u64 subsetnum(u64 x,int i) // The parameter i serves as a current line 
{
  	int       j; // current column 
  	u64       a; // modifies Pascal's coefficient a(n,j)
  	u64       p; // Pascal's coefficient 
  	u64	  w; // result (a bitmap)

  	// [Initialization]
  	w=0;
  	p=1;
  	a=0;
  	j=0;
	// The search for column j such that a(i,j) <= "x" < a(i,j+1)
  	// The Initialization of a(i,j)
  	while (a+p<=x)
  	{
    		a+=p;
    		j++;
    		p=(p*(i-j+1))/j;
  	}
  	p=(p*j)/i;
  	a=(a+p)/2;
  	--j;
  	x-=a;
  	--i;
 
  	while ( (j>=0) && (j<i) )
  	{ // The algorithm stops when we arrive on one of the sides. 
    	// The side j=0, in this case we complete with zéros
    	// The side i=j, we complete with 1 
    		if (x>=a) // Can we subtract a(i-1,j)
    		{
      			p=(p*(i-j))/i;
      			a=(a+p)/2; // Update of a
      			--i; // At the top only
    		}
    		else
    		{
      			w|=(1ll<<i); // remember that the number is 1
      			a=(a-(p*(i-j))/i)/2; // Update of a
      			p=(p*j)/i;  // Update of p
      			--i; // At the top 
     	 		--j; // At the left
    		}
    		x-=a; // Substracting what we could
  	}
  	// Fill last characters whith 0 or 1. 

  	if (j!=-1) 
  	{
    		w|=(1ll<<(i+1))-1ll;
  	}
  	return w;
}
  
// The function binome: return {n\choose p}
u64 binome(int n, int p)
{
    u64 r=1;
    int i;
    for (i=1;i<=p;i++,n--)
    {
        r=(r*n)/i;
    }
    return r;
}

#if (RANG==2)
// The group law C_f + C_e
// e is MODULO=EXPONENT of the group
// ------------------------
unsigned add(unsigned x,unsigned y)
{
	unsigned u,v;
	u=x%MODULO1;
	v=x/MODULO1;
	x=y%MODULO1;
	y=y/MODULO1;
	u+=x;
	if (u>=MODULO1) u-=MODULO1;
 	v+=y;
	if (v>=MODULO2) v-=MODULO2;
	return u+MODULO1*v;
}

// Display 
void affiche(unsigned x)
{
	printf("(%d,%d)",x%MODULO1,x/MODULO1);
}



// Opposite
unsigned opp(unsigned x)
{
	unsigned u,v;
	u=(x%MODULO1);
	if (u!=0) u=MODULO1-u;
	v=(x/MODULO1);
	if (v!=0) v=MODULO2-v;
	return u+MODULO1*v;
}
#endif


#if (RANG==3)

unsigned add(unsigned x, unsigned y)
{
	unsigned u,v,a,b;
	u=x%MODULO1;
	x/=MODULO1;
	a=y%MODULO1;
	y/=MODULO1;
	v=x%MODULO2;
	x/=MODULO2;
	b=y%MODULO2;
	y/=MODULO2;

	u+=a;
	if (u>=MODULO1) u-=MODULO1;
	v+=b;
	if (v>=MODULO2) v-=MODULO2;
	x+=y;
	if (x>=MODULO3) x-=MODULO3;
	return (v+x*MODULO2)*MODULO1+u;
}

unsigned opp(unsigned x)
{
	unsigned u,v;
	u=x%MODULO1;
	x/=MODULO1;
	v=x%MODULO2;
	x/=MODULO2;
	if (u!=0) u=MODULO1-u;
	if (v!=0) v=MODULO2-v;
	if (x!=0) x=MODULO3-x;
	return (v+x*MODULO2)*MODULO1+u;
}

void affiche(unsigned x)
{
	unsigned u,v;
	u=x%MODULO1;
	x/=MODULO1;
	v=x%MODULO2;
	x/=MODULO2;
	printf("(%u,%u,%u)",u,v,x);
}

#endif

#if (RANG==4)

unsigned add(unsigned x, unsigned y)
{
	unsigned u,v,w,a,b,c;
	u=x%MODULO1;
	x/=MODULO1;
	v=x%MODULO2;
	x/=MODULO2;
	w=x%MODULO3;
	x/=MODULO3;

	a=y%MODULO1;
	y/=MODULO1;
	b=y%MODULO2;
	y/=MODULO2;
	c=y%MODULO3;
	y/=MODULO3;

	u+=a;
	v+=b;
	w+=c;
	x+=y;
	if (u>=MODULO1) u-=MODULO1;
	if (v>=MODULO2) v-=MODULO2;
	if (w>=MODULO3) v-=MODULO3;
	if (x>=MODULO4) x-=MODULO4;
	return (((x*MODULO3+w)*MODULO2+v)*MODULO1+u);

}

unsigned opp(unsigned x)
{
	unsigned u,v,w;

	u=x%MODULO1;
	x/=MODULO1;
	v=x%MODULO2;
	x/=MODULO2;
	w=x%MODULO3;
	x/=MODULO3;
	if(u!=0) u=MODULO1-u;
	if(v!=0)v=MODULO2-v;
	if (w!=0) w=MODULO3-w;
	if (x!=0) x=MODULO4-x;
	return (((x*MODULO3+w)*MODULO2+v)*MODULO1+u);
	
}
#endif

//Calculate the opposite of the elements of a subset given by a bitmap. 
u64 completeoppsum(u64 w)
{
	u64 m;
	unsigned i,s;

	s=0;
	for (m=1,i=0;m<=(1ll<<CARDINAL);m<<=1,i++)
	{
		if ((m&w)!=0)
		{
			s=add(s,i);
		}
	}
	s=opp(s);
	m=1ll<<s; //
	if ((m&w)==0ll)
	{
		return m|w;
	}
	return 0ll;
}

// Calculate the sum of elements of a subset given by its bitmap w.
unsigned sum(u64 w)
{
	u64 m;
	unsigned s,i;
	for (s=0,i=0,m=1;m<(1ll<<CARDINAL);i++,m<<=1)
	{
		if ((m&w)!=0) s=add(s,i);
	}
	return s;
}

// The function w64: return the binary weight of x

int w64(u64 x)
{
	x=(x&0x5555555555555555ll)+((x>>1)&0x5555555555555555ll);
	x=(x&0x3333333333333333ll)+((x>>2)&0x3333333333333333ll);
	x=(x&0x0f0f0f0f0f0f0f0fll)+((x>>4)&0x0f0f0f0f0f0f0f0fll);
	x=(x&0x00ff00ff00ff00ffll)+((x>>8)&0x00ff00ff00ff00ffll);
	x=(x&0x0000ffff0000ffffll)+((x>>16)&0x0000ffff0000ffffll);
	return (x&0x00000000ffffffff)+((x>>32)&0x00000000ffffffffll);
}

// The function include: add an element in a set. 
void include(u64 x, u64*S)
{
	((u8*)S)[x>>3]|=(1<<(x&7));
}

// x is it in S?
int isin(u64 x, u64*S)
{
	return (((u8*)S)[x>>3]>>(x&7))&1;
}

// The function: clear_set: Makes the set S of cardinal n empty 
void clear_set(u64*S,u64 n)
{
	u64 w=(n+63)>>6;
	u64 i;
	for(i=0;i<w;i++) S[i]=0;
}

u64 cardinal(u64*S,u64 n)
{
	u64 w=(n+63)>>6;
	u64 c=0;
	u64 i;
	for (i=0;i<w;i++) c+=w64(S[i]);
	return c;
}

void affiche_set(u64*S,u64 n)
{
	u64 w=(n+63)>>6;
	u64 i;
	for (i=0;i<w;i++) 
	{
		printf("%llx",S[i]);
	}
	printf("\n");
}



u64 *setA,*setB;

// Binomial coefficient table. 
u64 b[CARDINAL];
u64 c[CARDINAL];


// b[n,p]= b[n*(n+1)/2 + p]
u64 bb[(CARDINAL+2)*(CARDINAL+1)/2];
u64 cc[(CARDINAL+2)*(CARDINAL+1)/2];

// Initialization of tables bb and cc
void initbbcc()
{
	int i,j,k;
	for (i=0;i<=CARDINAL;i++)
	{
		k=i*(i+1)/2;
		bb[k]=1;
		for (j=1;j<i;j++)
		{
			bb[k+j]=bb[k-i+j]+bb[k-i+j-1];
		}
		bb[k+i]=1;

	}
	for (i=0;i<=CARDINAL;i++)
	{
		k=i*(i+1)/2;
		cc[k]=1;
		for (j=1;j<i;j++)
		{
			cc[k+j]=cc[k-i+j]+cc[k-i+j-1];
		}
		cc[k+i]=cc[k+i-1]+1;

	}

}

// return the subset of a given number using the precalculated tables

u64 subsetnum1(u64 x,int i) // the i parameter indicates the current line
{
  	int       j; // the current column
  	u64       a; // modified Pascal's Coefficient a(n,j)
  	u64       p; // Pascal's triangle Coefficient
  	u64	  w; // result bitmap
	int k;
  	// initialization
  	w=0;
  	p=1;
  	a=0;
  	j=0;
  	// lookup for the column j such that a(i,j) <= "x" < a(i,j+1)
  	// et initialisation de a(i,j)
	k=i*(i+1)/2;
  	while (a+p<=x)
  	{
    		a+=p;
    		j++;
		p=bb[k+j];
  	}
  	--j;
	k-=i;
  	--i;
	a=cc[k+j];
	x-=a;
  	// nominale lookup
  	while ( (j>=0) && (j<i) )
  	{ // the algorithm stops whenever a side is reached
    	// let be the side j=0, then we complete with zeros
    	// let be the side i=j, then we complete with ones 
    		if (x>=a) // can we substruct a(i-1,j)
    		{
      			k-=i;
			--i; // to the top only
			a=cc[k+j];
    		}
    		else
    		{
      			w|=(1ll<<i); // memorize that the number is 1
			k-=i;
      			--i; // to the top
     	 		--j; // to the left
			a=cc[k+j];
    		}
    		x-=a; // we substruct with what we could
  	}
  	// filling the last characters with 0 or 1
  	if (j!=-1)
  	{
    		w|=(1ll<<(i+1))-1ll;
  	}
  	return w;
}
  
// return the number of a subsets using the precalculated tables
u64 numsubset1(u64 w,int n)
{
   	int i,j,k; // index of the course of the triangle
   	u64 a;   // pascal's triangle values
   	u64 x;  // cumul of all browsed summits 
   	// initialization
   	i=-1;
   	j=-1;
   	a=0;
   	x=0;
	k=0;
   	while (n>0)
   	{
      	// depends on the value of the number
     		i++; // move down anyways
		k+=i;
		switch (w&1)
      		{
        	case 0:
          	//	++i;           // move just down
          		break;
        	case 1:
          		//++i; 
			++j; // Move diagonally down right 
          		break;
      		}
      		// chift w to position the next number
      		w>>=1;
      		// cumulative value of the summit
		if (j==-1) a=0; else a=cc[k/*i*(i+1)/2*/+j];

	      	x+=a; //      		
		--n; // lookup counter
   	}
   	return x;
}



// calculate the list of next elements of a given number
// return the number of elements in the list
// the indexes are affected to *s
int suivants(u64 x,u64*s,int i) 
{
	int 	l; // result 
	int	j; // the current column
  	u64	a; // the modified pascal's triangle Coefficients a(n,j)
  	u64	p; // pascal's triangle Coefficients
	int k;
	u64 c;
  	// initialization
  	p=1;
  	a=0;
  	j=0;
	l=0;
  	// lookup for j column j such as a(i,j) <= "x" < a(i,j+1)
  	// and initialization of a(i,j)
	k=i*(i+1)/2;
  	while (a+p<=x)
  	{
    		a+=p;
    		j++;
		p=bb[k+j];
  	}
  	--j;
	k-=i;
  	--i;
	a=cc[k+j];
	c=+x;
	x-=a;
  	// nominal lookup
  	while ( (j>=0) && (j<i) )
  	{ // the algorithm stops whenever one side is reached
    	// let be the side j=0, then we complete with zeros
    	// let be the side i=j, then we complete with ones
    		if (x>=a) // can we substract a(i-1,j)
    		{// number 0, must memorize the next element
			c=s[l++]=c+bb[k+j+1];
      			k-=i;
			--i; // just at top
			a=cc[k+j];
    		}
    		else
    		{//  number 1, no next element, just memorize 
			// the value to be incremented
			c+=bb[k+j+1];
			k-=i;
      			--i; // at the top 
     	 		--j; // at the left 
			a=cc[k+j];
    		}
    		x-=a; // We substract what we could 
  	}
	if (j==-1)// If we stop at j==-1, il means that the bitmap ends with 0, and we have to add all the next elements. 
	{
		while (i>=0)
		{
			c=s[l++]=c+1;
			i--;
		}
	}
  	return l;
}
 

// Fill the set number k in set1 from date of set2. 
u64 brol(u64*set1,u64*set2,int k)
{
	u64 ii,pp,w,n;
	printf("p=%d b[p]=%llu\n",k,b[k]); fflush(stdout);
	clear_set(set1,b[k]);
	for (ii=0;ii<b[k-1];ii++)
	{
		if (isin(ii,set2)!=0)
		{
			w=subsetnum1(c[k-2]+ii,CARDINAL);
			for (pp=1;pp<(1ll<<CARDINAL);pp<<=1)
			{
				if ((pp&w)==0)
				{
					n=numsubset1(w|pp,CARDINAL);
					include(n-c[k-1],set1);

				}
			}
		}
	}
	pp=cardinal(set1,b[k]);
	printf("%llu %llu %llu\n",pp, b[k], b[k]-pp);
	fflush(stdout);
	return b[k]-pp;
}


// The procedure to obtain the successors  
u64 brol1(u64*set1,u64*set2,int k)
{
	u64 ii,pp;
	u64 s[CARDINAL];
	int l,i;

	printf("p=%d b[p]=%llu\n",k,b[k]); fflush(stdout);
	clear_set(set1,b[k]);
	for (ii=0;ii<b[k-1];ii++)
	{
		if (isin(ii,set2)!=0)
		{
			l=suivants(c[k-2]+ii,s,CARDINAL);
			
			for (i=0;i<l;i++)
			{
				include(s[i]-c[k-1],set1);
			}
		}
	}
	pp=cardinal(set1,b[k]);
	printf("%llu %llu %llu\n",pp, b[k], b[k]-pp);
	fflush(stdout);
	return b[k]-pp;
}


// The function verif: verifying that the sum of the elements of a subset is zero. 
// //------------------------------------------------------------------
void verif(u64 w)
{
	int s,i;
	s=0;
	u64 m;
	for (i=0,m=1;m<(1ll<<CARDINAL);m<<=1,i++)
	{
		if ((w&m)!=0) s=add(s,i);
	}
	if (s!=0) printf("!!! "); fflush(stdout);
}



// ************************* The algorithme description ************************* 

// Let : G  a finit abelian group ;
// N :  The cardinal of G ; 
// E : The  exponent of G ;
// g : Harboth’s constant.

//— The starting point: We consider all subsets of E elements (we have N choose E subsets of E elements).
//— We mark all of them whose sum’s zero.
//— Thus we obtain all N choose E subsets of E elements whose sum’s zero.
//— We browse all subsets of E elements. When one subset is marked, we mark
//all its successors immidiate for inclusion those how are obtained by adding any element.
// Thus we obtain all subsets of E + 1 elements which have a zero sum subset of E
//elements.
//— We repeat recursively with the subsets of E+1 elements marked until all the susets would be marked with a given cardinal.
//_ This cardinal represent the Harborth constant g(G).  
//— We note that just before reaching this constant, the subset no marked have no zero sum subset of E elements, so we can have ideas for the lower bound.
//— The subsets of G are represent by a bitmap.

 


int main()
{
	clock_t h;
	double dh;
	u64 no,m,ii,jj;
	u64 w;
	int n,p;
	printf("Hello somme nulle!\n");
	setA=malloc(3500000000);
	setB=malloc(3500000000);
	printf("setA size = %lu\n", malloc_usable_size(setA));
	printf("setB size = %lu\n", malloc_usable_size(setB));
	initbbcc();
	int k;

	// binominal coefficients table initialization
	n=CARDINAL;
	for(p=0;p<CARDINAL;p++) b[p]=binome(n,p);
	// calculate the cumulates
	c[0]=b[0];
	for(p=1;p<CARDINAL;p++) c[p]=c[p-1]+b[p];

	for (p=0;p<=2*EXPOSANT+1;p++)
	{
		printf("binome(%d,%d)=%llu cumuls=%llu\n",CARDINAL,p,b[p],c[p]);
	}

	no=c[EXPOSANT-2];
	m=c[EXPOSANT-1];
	printf("%llu ... %llu\n",no,m);

	
	// initialization for k=EXPONENT
	clear_set(setA,b[EXPOSANT]);
	jj=0;
	// browse of all sets of EXPONENT-1 elements
	for (ii=0;ii<b[EXPOSANT-1];ii++)
	{
		w=subsetnum1(ii+c[EXPOSANT-2],CARDINAL);

		w=completeoppsum(w);

		if (w!=0) jj++;
		if(w!=0)
		{

			include(numsubset1(w,CARDINAL)-c[EXPOSANT-1],setA);
		}
	}

	printf("%llu %llu %llu %llu\n",b[EXPOSANT-1],m-no,jj,cardinal(setA,b[EXPOSANT]));
	fflush(stdout);
	
	u64*s1,*s2,*st;
	s1=setB;
	s2=setA;
	u64 cs;
	k=EXPOSANT+1;
	do
	{	
		cs=brol1(s1,s2,k);
		st=s1;
		s1=s2;
		s2=st;
/* Listing of elements of interesting cases, in fact 
just before reaching this constant, the subset no marked have no zero sum subset of E elements, so we can have ideas for the lower bound.	
		if (k==10)
		{
			u64 ii;
			printf("sous-ensembles sans somme nulle:\n");
			for (ii=0;ii<b[k];ii++)
			{
				if (isin(ii,st)==0)
				{
					int l;
					u64 w;
					w=subsetnum1(ii+c[k-1],CARDINAL);
					for (l=0;l<CARDINAL;l++)
					{
						if (((1ll<<l)&w)!=0)
						{
							affiche(l);
							printf(" ");
						}
					}
					printf("\n");
					fflush(stdout);
				}
			}
		}
*/

		k++;
	}
	while ((cs!=0)&&(k<CARDINAL-1));

	h=clock();
	dh=((double)h)/CLOCKS_PER_SEC;
	printf("temps passé = %lf\n",dh);
	return 0;
}



