/* structure for zones data */

typedef struct region{
       float x1;
       float x2;
       float y1;
       float y2;
       float zdmax;
}REG;

typedef struct zones {
   int nz;
   REG *zlat;
   float zonemax;
}ZONE;

typedef struct adapt {
   float ct[4];
   float phii[4];
   float psl[3];
   float incp[3];
   float maxad;
}ADPT;
