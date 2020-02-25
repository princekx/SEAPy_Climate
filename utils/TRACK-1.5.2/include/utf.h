/* UTF defines */


#define  ARNG     3843.0         /* data range split into CNUM*CNUM-1 bins */
#define  CNUM     62
#define  MHD      15
#define  NFT      32
#define  RECL     2*NFT
#define  FTOL     1.0e-10
#define  CHMX     20
#define MAXCHR    500
#define  EOD      "***END OF DATA***"
#define CMISS     '?'
#define VMISS     MVAL

static char lkup[63] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz ";


char read_ascii( FILE * );
