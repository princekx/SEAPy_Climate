/* structure for boundary representation */

struct boundary_pt {
    union u_pt x;
    union u_pt y;
    union u_pt z;
    float val;
};


struct boundary_cntl {
    int fd;
    int nadd;
    int nmode;
    int ofill;
    int obex;
    int nexp;
};
