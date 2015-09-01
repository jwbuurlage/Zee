// header:
// int max_u
// int max_v
// int max_A
// int max_non_local_v
typedef a_header_t int[4];

typedef v_header_t int;

int main()
{
    a_header_t a_header;
    ebsp_get_header(0, &header);

    v_header_t v_header;
    ebsp_get_header(1, &v_header);

    float* u = malloc(max_u);
    float* v = malloc(max_u);
    float* vlocal = malloc(max_u);
    float* block = malloc(max_u);
    float* nonlocal_idxs = malloc(max_u);
    float* nonlocal_owners = malloc(max_u);

    // FIXME horizontal_strips
    for (v_strip < horizontal_strips) {
        ebsp_next_chunk(1);

        // FIXME vert_blocks
        for (int i = 0; i < vert_blocks; ++i) {
            ebsp_next_chunk(0);

            float* remote_idxs;
            float* remote_owners;

            for (;;) {
                bsp_hpget();
            }

            for (triplets) {
                // FMADD vs compressed
                u[i] += value * v[j];
            }

            // send u up
            // reset u, v
        }
    }


    return 0;
}
