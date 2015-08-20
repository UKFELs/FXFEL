Dear Jonathan,

Here is a version of the MASP software. MASP stands for
"MAcroparticle Splitting Program", and its name is
provisional.

There are also two pdf papers. The Specification explains
how to use the program, and the other paper explains
how it works. There is some overlap between the two
papers, but not much.

To run the system, load all the .c and .h modules into
a working directory, together with  TOPHAT,BLANCMANGE
and TRIFLE. You will also need a data set. I supply
three:

electrons.txt, with about 38000 records 
          (use DATA 0 1 2 3 4 5 7)
german.txt  with about 11000 records  
          (use DATA 0 3 1 4 2 5 11)
reduced.txt with about 1900 records
This last is a random sample of electrons.txt
          (use DATA 0 1 2 3 4 5 7)
      
To compile and run the check program, I use

gcc maspcheck.c -o zz -lm
./zz (parameter file)

To compile and run the main program I use
mpicc masp.c -o pa -lm
mpirun -np 12 ./pa (parameter file)

If you need more cores, then on Archie I would
use a job script, as quoted in the specification.
I gather your environment is different, but then
you would know what to do.

This is the first time I have let the program 
out of my control, and I expect numerous comments,
bug alerts, etc. If you let me know what they are
I will attend to them promptly.

Please do not make significant changes to the code
without consulting me. 

You can contact me in several ways:

Land line : 0141 956 6746
Mobile:	    07913562860
Email:	    andrew@crm.scotnet.co.uk
SKYPE :     andrew.colin7  (I think)

