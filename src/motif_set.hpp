#ifndef __MOTIF_SET__
#define __MOTIF_SET__

#include "motif.hpp"

#define MOTIF_SET_SIZE 1024

class MotifSet
{
        typedef tuple<const Motif*, uint64, uint64> elem_desc_t; //motif_array, start_pos, end_pos
        typedef pair<Motif, uint32> heap_elem_t; //motif val, desc_id
        elem_desc_t data_desc[MOTIF_SET_SIZE];
        heap_elem_t data[MOTIF_SET_SIZE];
        uint32 pos;
        uint32 desc_pos;

        inline void update_heap()
        {
                uint32 desc_id = data[1].second;
                Motif motif;
                if (++get<1>(data_desc[desc_id]) < get<2>(data_desc[desc_id]))
                {
                        motif.set(get<0>(data_desc[desc_id])[get<1>(data_desc[desc_id])]);
                }
                else
                {
                        motif.set(data[--pos].first);
                        desc_id = data[pos].second;
                }

                uint32 parent, less;
                parent = less = 1;
                while (true)
                {
                        if (parent * 2 >= pos)
                                break;
                        if (parent * 2 + 1 >= pos)
                                less = parent * 2;
                        else if (data[parent * 2].first < data[parent * 2 + 1].first)
                                less = parent * 2;
                        else
                                less = parent * 2 + 1;
                        if (data[less].first < motif)
                        {
                                data[parent] = data[less];
                                parent = less;
                        }			
                        else
                                break;
                }
                data[parent] = make_pair(motif, desc_id);
        }

        public:
        MotifSet()
        {
                pos = 1;
                desc_pos = 0;
        }
        inline void init_add(const Motif* buffer, uint64 start_pos, uint64 end_pos)
        {
                data_desc[desc_pos] = make_tuple(buffer, start_pos, end_pos);
                data[pos].first.set(buffer[start_pos]);
                data[pos].second = desc_pos;
                uint32 child_pos = pos++;

                while (child_pos > 1 && data[child_pos].first < data[child_pos / 2].first)
                {
                        swap(data[child_pos], data[child_pos / 2]);
                        child_pos /= 2;
                }
                ++desc_pos;
        }
        inline void clear()
        {
                pos = 1;
                desc_pos = 0;
        }

        inline bool get_min(Motif& motif)
        {
                if (pos <= 1)
                        return false;

                motif = data[1].first;
                update_heap();


                return true;
        }
};



#endif // __MOTIF_SET__
