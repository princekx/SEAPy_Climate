/* header file for the image structure */

struct eq_class {
    struct eq_class *farther;
    int name;
};


struct image {
     int pval;
     struct image *mother;
     struct image *nw;
     struct image *ne;
     struct image *sw;
     struct image *se;
     struct eq_class *label;
     struct eq_class *stor;
};
