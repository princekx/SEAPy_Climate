/* structure definitions for tracks */

struct track_points{
        int frame_id;
        int object_id;
        int feature_id;
        int nmpt;
};


struct track_ind{
        int num_point_track;
        struct track_points *tp;
};
