#ifndef MY_ISNAN_H
#define MY_ISNAN_H

#ifdef BIG_ENDIAN_ON
typedef union 
{
  double value;
  struct 
  {
    unsigned int sign : 1;
    unsigned int exponent: 11;
    unsigned int fraction0:4;
    unsigned int fraction1:16;
    unsigned int fraction2:16;
    unsigned int fraction3:16;
    
  } number;
  struct 
  {
    unsigned int sign : 1;
    unsigned int exponent: 11;
    unsigned int quiet:1;
    unsigned int function0:3;
    unsigned int function1:16;
    unsigned int function2:16;
    unsigned int function3:16;
  } nan;
  struct 
  {
    unsigned long msw;
    unsigned long lsw;
  } parts;
    long aslong[2];
} ieee_double_type;
#endif

#ifdef LITTLE_ENDIAN_ON
typedef union 
{
  double value;
  struct 
  {
    unsigned int fraction3:16;
    unsigned int fraction2:16;
    unsigned int fraction1:16;
    unsigned int fraction0: 4;
    unsigned int exponent :11;
    unsigned int sign     : 1;
  } number;
  struct 
  {
    unsigned int function3:16;
    unsigned int function2:16;
    unsigned int function1:16;
    unsigned int function0:3;
    unsigned int quiet:1;
    unsigned int exponent: 11;
    unsigned int sign : 1;
  } nan;
  struct 
  {
    unsigned long lsw;
    unsigned long msw;
  } parts;
  long aslong[2];
} ieee_double_type;
#endif




inline bool myisnan (double x)
{
  ieee_double_type val;
  val.value = x;
//   fprintf (stderr, "sign      = %x\n", val.number.sign    );
//   fprintf (stderr, "exponent  = %x\n", val.number.exponent);
//   fprintf (stderr, "fraction0 = %x\n", val.number.fraction0);
//   fprintf (stderr, "fraction1 = %x\n", val.number.fraction1);
//   fprintf (stderr, "fraction2 = %x\n", val.number.fraction2);
//   fprintf (stderr, "fraction3 = %x\n", val.number.fraction3);
  if (((val.number.exponent & 0x7ff) == 0x7ff) &&
      ( (val.number.fraction0 !=0)  || 
	(val.number.fraction1 != 0) || 
	(val.number.fraction2 != 0)))
    return (true);
  else
    return false;
}


// #include <iostream>
// using namespace std;
	

// main ()
// {
//   double x = 3.0;
//   cerr << "x = " << x << endl;
//   cerr << "myisnan = " <<  (myisnan(x) ? "true" : "false")  << endl;
//   x = sqrt (-1.0);
//   cerr << "x = " << x << endl;
//   cerr << "myisnan = " <<  (myisnan(x) ? "true" : "false")  << endl;
// }
  
#endif
