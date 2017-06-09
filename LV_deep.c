#include "LV_deep.h"

// #define DEBUG
#ifdef DEBUG
#define debug_print(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
#else
#define debug_print(fmt, ...) do {} while (0)
#endif

#define LOGA 0

void init_node_space(DI* di, size_t s_l, int n_l)
{
	(*di).s_l = s_l;
    (*di).n_l = n_l;
    (*di).the_deep = -1;
    (*di).deep_s = (char *)malloc(sizeof(char) * (s_l + 10));
    (*di).all_node = (deep_node *)malloc(sizeof(deep_node) * (n_l));
    return;
}

void free_node_space(DI* di)
{
	free((*di).all_node);
	free((*di).deep_s);
	// free(di);
    return;
}



void mark_node_id(DI *di, int node_index, int deep, size_t begin, size_t len)
{
	(*di).all_node[node_index].id = node_index;
	(*di).all_node[node_index].deep = deep;
	(*di).all_node[node_index].begin = begin;
	(*di).all_node[node_index].len = len;
	(*di).the_deep = MX((*di).the_deep, deep+1);
	return;
}


void log_DI(DI *di)
{
	printf(" %d ", (*di).n_l);
	if (LOGA == 0) 
	return;
	// printf("s_l: %lu n_l: %d the_deep: %d\n", (*di).s_l, (*di).n_l, (*di).the_deep);
	// for (size_t i=0; i<(*di).s_l; ++i)
	// {
	// 	print_bit4( (*di).deep_s[i] );
	// }
	printf("\n");
	if ((*di).n_l >= 1)
	for (size_t i=0, k; i<(*di).n_l; ++i)
	{
		printf("deep: %d node: %d: ", (*di).all_node[i].deep, (*di).all_node[i].id);
		for (size_t t, j = 0; j<(*di).all_node[i].len; ++j)
		{
			t = (*di).all_node[i].begin + j;
			print_bit4( (*di).deep_s[t] );
		}
		printf("\n");
	}
}


//parse the target string into deep node and string
//it have a tree struct and a for ref deep string
void parse_target_string(DI *di, size_t *ll, size_t *rr, ref_s *rs, var_s *vs)
{
	size_t l = *ll, r = *rr;
	if (r<=l) return;
	l = ( l < rs->ref_s_l ) ? l : rs->ref_s_l;
	l = l/(vs->s_m)*(vs->s_m);
	r = ( r < rs->ref_s_l ) ? r : rs->ref_s_l;
	// *ll = l;
	*rr = r;

	//tmp
	size_t tar_n, tar_u, tmp;

	//index for retrive snp and mnp
	size_t s_cnt = vs->bv1_s[l / vs->s_m], m_cnt = vs->bv2_s[l / vs->s_m], begin_s_cnt, begin_m_cnt;

	//read the string before l
	for (size_t k = l; k < *ll; ++k)
	{
		if (get_mbit(vs->bit_v1, k))
		{
			++s_cnt;
		}
		else if (get_mbit(vs->bit_v2, k))
		{
			get_HAbit(&m_cnt, vs->mnp, &tar_n, 0);
			for (int i=0; i<tar_n; ++i)
			{
				get_HAbit(&m_cnt, vs->mnp, &tar_u, 1);

				
				for (int j=0; j<tar_u; ++j)
				{
					get_Abit(&m_cnt, vs->mnp, &tmp, 2);
				}

			}
		}
	}

	begin_s_cnt = s_cnt;
	begin_m_cnt = m_cnt;

	//for init space

	size_t s_l=0, n_l = 0;

	//			  o
	// 			/ 	\
	//		o - 	  - o
	//			\   /
	//            o
	//  the number of node, node len is 4
	for (size_t k = *ll; k < r; ++k)
	{
		if (get_mbit(vs->bit_v2, k))
		{
			// get_Abit(&m_cnt, vs->mnp, &tar_n, alt_n_l);
			get_HAbit(&m_cnt, vs->mnp, &tar_n, 0);
			n_l += (tar_n + 1);
			for (int i=0; i<tar_n; ++i)
			{
				// get_Abit(&m_cnt, vs->mnp, &tar_u, alt_u_l);
				get_HAbit(&m_cnt, vs->mnp, &tar_u, 1);

				s_l += tar_u;
				if (i==0)
				{
					k += MX((tar_u-1), 0);
				}
				for (int j=0; j<tar_u; ++j)
				{
					get_Abit(&m_cnt, vs->mnp, &tmp, 2);
				}
			}
		}
		else 
		{
			s_l += 1;
		}
	}
	init_node_space(di, s_l, n_l+1);

    s_cnt = begin_s_cnt;
	m_cnt = begin_m_cnt;


	size_t last_b = 0, d_s_i=0;
	size_t node_index = 0;
	int deep = 0;

	for (size_t k = *ll; k < r; ++k)
	{
		if (get_mbit(vs->bit_v2, k))
		{
			get_HAbit(&m_cnt, vs->mnp, &tar_n, 0);

			mark_node_id(di, node_index++, deep++, last_b, d_s_i - last_b);
			last_b = d_s_i;

			// printf("mnp pos: %lu #%lu\n", k+1, tar_n);
			for (int i=0; i<tar_n; ++i)
			{
				// get_Abit(&m_cnt, vs->mnp, &tar_u, alt_u_l);
				get_HAbit(&m_cnt, vs->mnp, &tar_u, 1);
				if (i==0)
				{
					k += MX(0, tar_u-1);
				}
				for (int j=0; j<tar_u; ++j)
				{
					get_Abit(&m_cnt, vs->mnp, &tmp, 2);
					// fprintf(out, "%lu ", tmp);
					(*di).deep_s[d_s_i++] = map_bit4[map_c[tmp]];
				}
				mark_node_id(di, node_index++, deep, last_b, tar_u);
				last_b = d_s_i;
			}
			++deep;
		}
		else if (get_mbit(vs->bit_v1, k))
		{
			// printf("*");
			(*di).deep_s[d_s_i++] = get_snp(s_cnt, vs->snp) | map_bit4[get_ref_c(rs, k)];
			++s_cnt;
		}
		else 
		{
			(*di).deep_s[d_s_i++] = map_bit4[get_ref_c(rs, k)];
		}

	}

	mark_node_id(di, node_index++, deep++, last_b, d_s_i - last_b);
	last_b = d_s_i;

	for (size_t k = last_b; k < (last_b+9); ++k)
	{
		(*di).deep_s[k++] = 0;
	}


	log_DI(di);
    return;
}


//print a bit 4 to terminal.ie 0011 should print [AC]
void print_bit4(char o)
{
	char  z=0;
	char tmp[10];
	memset(tmp, 0, sizeof(tmp));
	if (o & 1) tmp[z++] ='A';
	if (o & 2) tmp[z++] ='C';
	if (o & 4) tmp[z++] ='G';
	if (o & 8) tmp[z++] ='T';
	if (z>1)	printf("[%s]", tmp);
	else printf("%s", tmp);
}

//print a bit 4 to terminal.ie 0011 should print [AC]
void print_bit4_(char o)
{
	char  z=0;
	char tmp[10];
	memset(tmp, 0, sizeof(tmp));
	if (o & 1) tmp[z++] ='A';
	if (o & 2) tmp[z++] ='C';
	if (o & 4) tmp[z++] ='G';
	if (o & 8) tmp[z++] ='T';
	if (z>1)	printf("*");
	else printf("%s", tmp);
}

//int lv space should large enough
void init_lv_space(LV_ENTITY *lv, int n_l, int max_error, int lp, int deep)
{
	int ***LV;
	LV = (int ***)malloc(sizeof(int **) * (n_l));
	for(size_t b=0; b<n_l; ++b)
	{
		LV[b] = (int **)malloc(sizeof(int *) * (max_error + 2));
		for(size_t e=0; e<max_error+2; ++e)
		{
			LV[b][e] = (int *)malloc(sizeof(int) * ((lp+1)*2 + 1));
		}
	}

	int **Far_reach;
	Far_reach = (int **)malloc(sizeof(int *) * (n_l));
	for(size_t b=0; b<(n_l); ++b)
	{
		Far_reach[b] = (int *)malloc(sizeof(int) * ((lp+1)*2 + 1));
	}

	int **mark_bit;
	mark_bit = (int **)malloc(sizeof(int *) * (deep));
	for(size_t dp=0; dp<(deep); ++dp)
	{
		mark_bit[dp] = (int *)malloc(sizeof(int ) * (lp + 1));
	}

	lv->n_l = n_l;
	lv->max_error = max_error;
	lv->lp = lp;
	lv->deep = deep;

	lv->LV = LV;
	lv->Far_reach = Far_reach;
	lv->mark_bit = mark_bit;
	return;
}

void free_lv_space(LV_ENTITY *lv)
{
	for(size_t b=0; b<(lv->n_l); ++b)
	{
		for(size_t e=0; e<(lv->max_error)+2; ++e)
			free((lv->LV)[b][e]);
		free((lv->LV)[b]);
	}
	free(lv->LV);
	for(size_t b=0; b<(lv->n_l); ++b)
		free((lv->Far_reach)[b]);
	free(lv->Far_reach);

	for(size_t b=0; b<(lv->deep); ++b)
		free((lv->mark_bit)[b]);
	free(lv->mark_bit);

	// free(di);
    return;
}



//init lv all mark_no_vis
int init_lv(LV_ENTITY *lv, int n_l, int max_error, int lp, int deep, int mark_no_vis)
{
	if (n_l > (lv->n_l) || max_error > (lv->max_error) || lp > (lv->lp))
	{
		printf("init lv error, space require out of range:\n");
		printf("number of block max : %d , require %d\n", (lv->n_l), n_l);
		printf("max error max : %d , require %d\n", (lv->max_error), max_error);
		printf("len of pattern max : %d , require %d\n", (lv->lp), lp);
		return -1;
	}

	for(size_t b=0; b<n_l; ++b)
	{
		for(size_t e=0; e<max_error+2; ++e)
		{
			// ek = e + 1;
			// dk = d + lp + 1;
			for(size_t d=0; d<=(lp+1)*2; ++d)
			{
				(lv->LV)[b][e][d] = mark_no_vis;
			} 
		}
	} 

	// log
	// for(size_t b=0; b<n_l; ++b)
	// {
	// 	for(size_t e=0; e<=max_error+1; ++e)
	// 	{
	// 		for(size_t d=0; d<=lp*2; ++d)
	// 		{
	// 			printf("%d", (lv->LV)[b][e][d]);
	// 		}
	// 		printf("\n");
	// 	}
	// 	printf("\n");
	// }


	for(size_t b=0; b<n_l; ++b)
	{
		for(size_t d=0; d<=lp*2; ++d)
		{
			(lv->Far_reach)[b][d] = 0;
		}
	}

	return 0;
}


//deep node string and pattern should parsed
void run_deep_lv(LV_ENTITY *lv, DI *di, char *pattern, size_t pattern_len, int max_error)
{
	int mark_no_vis = -1;

	//lp should inside [INT] range
	if (pattern_len > INT_MAX)
	{
		fprintf(stderr, "len of pattern [%lu] out of range try small len\n", pattern_len);
		return;
	}
	int lp = (int)pattern_len;

	if (init_lv(lv, (di->n_l), max_error, lp, (di->the_deep), mark_no_vis) == -1)
	{
		fprintf(stderr, "run lv error\n");
		return;
	}
	

    int find_d = 0;
	int global_reach = -2;

	int the_b = 0, the_e = 0, the_d = 0;
	int the_ek = the_e + 1, the_dk = the_d+lp+1;

	// init the LV[0][0][0] = 0
	// (lv->LV)[the_b][the_ek][the_dk] = 0;

	//tmp var 
	deep_node *bn;
	int ek, dk, best, last_best;
	//len of text , text begin in deep code
	int lt;
	size_t tbegin;

	char *p, *pb, *t, *pend;
    int extend_len = 0;

    int the_deep = 0, base_deep, cur_deep, action, next_action;

	for (int e = 0; e<=max_error && find_d==0; ++e)
	{
		//mark_bit mean if lv == no vis 
		//mark_bit = 1 lv can succeed from d+1 or == -d then extend 0 error
		//mark_bit = 2 lv can succeed then extend 1 error
		for (int dp = 0; dp <= the_deep; ++dp)
		{
			for (int i=0; i<lp; ++i)
			{
				if (dp==0)
					if (i<=e) (lv->mark_bit)[dp][i] = e-i+1;
					else (lv->mark_bit)[dp][i] = 0;
				else
					(lv->mark_bit)[dp][i] = 0;
			}	
		}

		base_deep = -1;

		debug_print("\ne: %d\n", e);
		if ((di->n_l) > 0)
		// for (int b = 0; b<(di->n_l) && ((di->all_node[b]).deep <= the_deep); ++b)
		for (int b = 0; b<(di->n_l) && ((di->all_node[b]).deep <= the_deep) && find_d==0; ++b)
		{
			// if ((di->all_node[b]).deep == cur_deep) continue;
			bn = &(di->all_node[b]);
			lt = (int)MN(bn->len, INT_MAX);
			tbegin = bn->begin;
			cur_deep = bn->deep;

			debug_print("$ deep: %d the node : %d len: %d \n", cur_deep, b, lt);
			// printf("$ deep: %d the node : %d len: %d \n", cur_deep, b, lt);

			for (size_t i=0; i<lp; ++i)	debug_print((i==lp-1) ? "%d\n":"%d-", (lv->mark_bit)[cur_deep][i]);

			//if the block have no element,just succeed the mark value
			// if (lt == 0 && cur_deep+1 < (di->the_deep))
   //          {
   //          	if (the_deep < cur_deep+1)
   //          	{
   //          		the_deep = cur_deep+1;
   //          		for (size_t i=0; i<lp; ++i)	(lv->mark_bit)[cur_deep+1][i] = 0;
   //          		debug_print(" ** undate deep: cur %d the_deep %d\n", cur_deep, the_deep);
   //          	}
   //          	for (size_t i=0; i<lp; ++i)
   //          		(lv->mark_bit)[cur_deep+1][i] = MX((lv->mark_bit)[cur_deep+1][i], (lv->mark_bit)[cur_deep][i]);
   //          	// if ((lv->mark_bit)[cur_deep+1][best] > 1)
   //          	// {
   //          	// 	if (best+1 > 0 && best+1<lp) (lv->mark_bit)[cur_deep+1][best+1] = MX((lv->mark_bit)[cur_deep+1][best+1], 1);
   //          	// 	if (best-1 > 0 && best-1<lp) (lv->mark_bit)[cur_deep+1][best-1] = MX((lv->mark_bit)[cur_deep+1][best-1], 1);
   //          	// }
			// 	for (size_t i=0; i<lp; ++i)	debug_print((i==lp-1) ? "%d\n":"%d-", (lv->mark_bit)[cur_deep+1][i]);
   //          }
			

			for (int d = -lp+1, ek; d<(int)(MN( e+1, MN(lt+1, lp))); ++d)
			{
				ek = e + 1;
				dk = d + lp + 1;
				
				// if ((lv->Far_reach)[b][dk]==1) continue;
				if ((lv->Far_reach)[cur_deep][lt-d]==1) continue;

				best = MX((lv->LV)[b][ek-1][dk], (lv->LV)[b][ek][dk]);
				// debug_print("d: %d : best %d\n", d, best);

				// best =(lv->LV)[b][e ][dk];
				if (best != mark_no_vis) action = 2;
				else if (d <= 0) 
				{
					action = (lv->mark_bit)[cur_deep][-d];
					if (action > 0) 
					{
						// printf("-- %d %d %d\n", (lv->LV)[b][ek-1][dk+1], (lv->LV)[b][ek-1][dk-1], best);
						if ((lv->LV)[b][ek-1][dk+1] != mark_no_vis)
							best = MX(best, (lv->LV)[b][ek-1][dk+1] + 1);
						if ((lv->LV)[b][ek-1][dk-1] != mark_no_vis)
							best = MX(best, (lv->LV)[b][ek-1][dk-1]);
						best = MX(-d, best);
					}
				}
				else if ((lv->LV)[b][ek-1][dk-1] != mark_no_vis) action = 2;
				else action = 0;


				if (best == mark_no_vis && action == 0) continue;

				debug_print("action: %d\n", action);
				debug_print("d: %d : best %d\n", d, best);

				next_action = 0;
				
				//error 0 extend
				//only can extend from left side
				if (action >= 1 && (lv->LV)[b][ek-1][dk] == mark_no_vis && d<=0)
				{
					pb = pattern + best;
					p = pattern + best;
		    		t = di->deep_s + tbegin + best + d;
		    		last_best = best;

		    		extend_len = MN(lt-MX(0, d), lp-MX(0, -d));
		    		pend = pattern + MX(0, -d) + extend_len;
		    		// debug_print("extend_len, %d\n", extend_len);
					// printf("1t: ");
					// print_bit4((*p));
					// printf("-");
					// print_bit4((*t));
					if ( p < pend && ((*p) & (*t)))
					{
						while (p < pend) 
						{
			                uint64_t x = (~ *((uint64_t*) t) ) & *((uint64_t*) p);
		    		        if (x) {
		    		            unsigned long zeroes;
		    		            CountTrailingZeroes(x, zeroes);
		    		            zeroes >>= 3;
								debug_print("* %lu: ", zeroes);
		    		            best = MN((int)(p - pb) + (int)zeroes + last_best, MX(0, -d) + extend_len);
		    		            
		    		            break;
		    		        }
		    		        p += 8;
		    		        t += 8;
		    		        if (p >= pend) 
		    		        {
		    		        	debug_print("@ %d: ", extend_len);
		                        best = MX(0, -d) + extend_len;
		                        break;
		                    }
		    		    }
					}
					if (best != last_best)
					{
						(lv->LV)[b][ek][dk] = best;
					}
					// else
					// {
					// 	if (d <= 0 && best == -d)
					// 		best = mark_no_vis;
					// }

					if (best + d >= lt)
					{
						next_action = action;
					}
					
				}
				debug_print("--d: %d : best %d\n", d, best);

	    		//error 1 extend
	    		if (action >= 2 && next_action == 0)
	    		{
	    			if (best != mark_no_vis) 
	    				best = best + 1;

					if ((lv->LV)[b][ek-1][dk] != mark_no_vis)
						best = MX(best, (lv->LV)[b][ek-1][dk] + 1);

					if ((lv->LV)[b][ek-1][dk+1] != mark_no_vis)
						best = MX(best, (lv->LV)[b][ek-1][dk+1] + 1);

					if ((lv->LV)[b][ek-1][dk-1] != mark_no_vis)
						best = MX(best, (lv->LV)[b][ek-1][dk-1]);


					if (best == (lv->LV)[b][ek-1][dk] + 1 && best + d > lt)
					{
						next_action = 2;
					}
					if (next_action == 0 && best == (lv->LV)[b][ek-1][dk+1] + 1 && best + d > lt)
					{
						next_action = 1;
					}
					if (next_action == 0 && best == (lv->LV)[b][ek-1][dk-1] && best + d > lt)
					{
						next_action = 1;
					}
					debug_print("d: %d : best %d\n", d, best);

					pb = pattern + best;
					p = pattern + best;
		    		t = di->deep_s + tbegin + best + d;
		    		last_best = best;

		    		extend_len = MN(lt-MX(0, d), lp-MX(0, -d));
		    		// debug_print("extend_len, %d\n", extend_len);
		    		pend = pattern + MX(0, -d) + extend_len;
					
					// printf("t: ");
					// print_bit4((*p));
					// printf("-");
					// print_bit4((*t));

					if ( p < pend && ((*p) & (*t)))
					{
						// printf("try match");
						// print_bit4((*p));
						// printf("-");
						// print_bit4((*t));
						// printf("\n");
						while (p < pend) 
						{
			                uint64_t x = (~ *((uint64_t*) t) ) & *((uint64_t*) p);
		    		        if (x) {
		    		            unsigned long zeroes;
		    		            CountTrailingZeroes(x, zeroes);
		    		            zeroes >>= 3;
		    		            best = MN((int)(p - pb) + (int)zeroes + last_best, MX(0, -d) + extend_len);
		    		            
		    		            break;
		    		        }
		    		        p += 8;
		    		        t += 8;
		    		        if (p >= pend) 
		    		        {
		                        best = MX(0, -d) + extend_len;
		                        break;
		                    }
		    		    }
					}
					if (best != last_best)
					{
						(lv->LV)[b][ek][dk] = best;
					}

					if (next_action == 0 && best + d == lt)
					{
						next_action = 1;
					}

					if (next_action == 0 && best + d > lt)
					{
						next_action = 2;
					}
	    			
	    		}

				if (best == mark_no_vis) continue;

				best = MN(best, lt - d);
				
				debug_print("d: %d : best %d\n", d, best);


				// if (best <= (lv->LV)[b][ek-1][dk])
				if (best+d >= lt)
				{
					// (lv->Far_reach)[b][dk] = 1;
					(lv->Far_reach)[cur_deep][lt-d] = 1;
					debug_print("Far reach at b: %d deep: %d d: %d lv: %d\n", b, cur_deep, d, best);
				}


	            // if (best > global_reach || (best == global_reach && d > the_d) )
	            if (best > global_reach)
	            {
	            	global_reach = best;
	            	the_b = b;
	                the_e = e;
	                the_d = d;
	            }

	            (lv->LV)[b][ek][dk] = best;
				debug_print("-> b: %d e:%d d:%d lV:%d\n", b, e, d, best);

	            if (best + d >= lt && cur_deep+1 < (di->the_deep))
	            {
	            	if (the_deep < cur_deep+1)
	            	{
	            		the_deep = cur_deep+1;
	            		for (size_t i=0; i<lp; ++i)	(lv->mark_bit)[cur_deep+1][i] = 0;
	            		debug_print(" ** undate deep: cur %d the_deep %d\n", cur_deep, the_deep);
	            	}
	            	(lv->mark_bit)[cur_deep+1][best] = MX((lv->mark_bit)[cur_deep+1][best], next_action);
	            	// if ((lv->mark_bit)[cur_deep+1][best] > 1)
	            	// {
	            	// 	if (best+1 > 0 && best+1<lp) (lv->mark_bit)[cur_deep+1][best+1] = MX((lv->mark_bit)[cur_deep+1][best+1], 1);
	            	// 	if (best-1 > 0 && best-1<lp) (lv->mark_bit)[cur_deep+1][best-1] = MX((lv->mark_bit)[cur_deep+1][best-1], 1);
	            	// }
					for (size_t i=0; i<lp; ++i)	debug_print((i==lp-1) ? "%d\n":"%d-", (lv->mark_bit)[cur_deep+1][i]);
	            }

				
				if (best == lp)
	            {
	            	the_b = b;
	                the_e = e;
	                the_d = d;
	                if (find_d==1) break;
	                find_d = 1;
	                // printf(" *** find lv: in error: %d b: %d deep: %d d: %d\n", e, b, cur_deep, d);
	                // break;
	            }
			}
		}
		
	}

	if (LOGA)	printf("global_reach: e: %d lv:%d b: %d d: %d \n\n", the_e, global_reach, the_b, the_d);
	else printf("*%d ", the_e);
	// get_Cigar(lv, di, pattern, lp, max_error, \
	// 	global_reach,  the_b,  the_e,  the_d);
}


// void get_Cigar(LV_ENTITY *lv, DI *di, char *pattern, int lp, int max_error, \
// 		int global_reach, int the_b, int the_e, int the_d)
// {
// 	Cigar *d_path;
// 	d_path = (Cigar *)malloc(sizeof(Cigar) * (lp + the_e));

// 	int *deep_path;
// 	deep_path = (int *)malloc(sizeof(int) * (lp + the_e));

// 	int e=the_e, d=the_d, b=the_b, cur_deep;
// 	int dk, ek, the_n = global_reach-1, reach=-2;
// 	int prev_reach = 0;

// 	int d_index = 0;

// 	if ((di->n_l) == 0) 
// 	{
// 		size_t i = 0;
// 		while (i < max_error && i<lp)
// 		{
// 			printf("_");
// 			++i;
// 		}
// 		printf("\n");
// 		i = 0;
// 		while (i < max_error && i<lp)
// 		{
// 			int o = pattern[i];
//             print_bit4_(o);
//             ++i;
// 		}
// 		printf("\n");
// 		return;
// 	}
// 	deep_node *bn;
// 	size_t tbegin;
// 	bn = &(di->all_node[b]);
// 	tbegin = bn->begin;
// 	cur_deep = bn->deep;
// 	int lt = (int)MN(bn->len, INT_MAX);

// 	the_n = global_reach-1;
// 	while (e >= 0)
// 	{
// 		ek = e + 1;
// 		dk = d + lp + 1;
// 		// printf("e: %d b: %d d: %d n: %d \n", e, b, d, the_n);
				
// 		if ((d <= 0 && the_n+1 > -d) || (d>0 && the_n > 0))
// 		{
//             // printf("%d %d %d %d %d | ", e, d, (lv->LV)[b][ek][dk-1], (lv->LV)[b][ek][dk], (lv->LV)[b][ek][dk+1]);
// 			prev_reach = MX((lv->LV)[b][e][dk-1], MX((lv->LV)[b][e][dk+1]+1, (lv->LV)[b][e][dk]+1));
// 			// printf("prev_reach: %d \n", prev_reach-1);
			
// 			// prev_reach = MX((lv->LV)[b][e][dk-1], MX((lv->LV)[b][e][dk+1]+1, (lv->LV)[b][e][dk]+1));
// 			// if (prev_reach == (lv->LV)[b][e][dk]+1) prev_reach -= 1;
// 			if ((lv->LV)[b][e][dk]+1 <= the_n+1)
// 			{
// 				// printf("try match\n");
// 				int nn = the_n;
// 				int dd = d_index;
// 				while ((lv->LV)[b][e][dk]+1 <= the_n + 1 && the_n >= 0 && the_n+d >= 0 && \
// 					(pattern[the_n] & *(di->deep_s + tbegin + the_n + d)))
// 				{
					
// 					deep_path[d_index] = b;
// 					d_path[d_index++] = mat;
// 		            // printf ("mat | ");
// 		            // printf ("n: %d d: %d b: %d ", the_n, d, b);
// 		            // print_bit4(pattern[the_n]);
// 		            // printf ("| ");

// 		            --the_n;
// 				}
// 				if (the_n > (lv->LV)[b][e][dk])
// 				{
// 					the_n = nn;
// 					d_index = dd;
// 					if (prev_reach <= the_n + 1)
// 					{
// 						// printf("try match\n");
// 						// print_bit4(pattern[the_n]);
// 						// printf("-");
// 						// print_bit4(*(di->deep_s + tbegin + the_n + d));
// 						// printf("\n");

// 						while (prev_reach <= the_n + 1 && the_n >= 0 && the_n+d >= 0 && \
// 							(pattern[the_n] & *(di->deep_s + tbegin + the_n + d)))
// 						{
							
// 							deep_path[d_index] = b;
// 							d_path[d_index++] = mat;
// 				            // printf ("mat | ");
// 				            // printf ("n: %d d: %d b: %d ", the_n, d, b);
// 				            // print_bit4(pattern[the_n]);
// 				            // printf ("| ");

// 				            --the_n;
// 						}
// 					}
// 					// printf(" %d %d %d %d\n", the_n, nn, d_index-dd, prev_reach);
// 				}

// 			}
// 			// if (prev_reach <= the_n + 1)
// 			// {
// 			// 	printf("try match ");
// 			// 	print_bit4(pattern[the_n]);
// 			// 	printf("-");
// 			// 	print_bit4(*(di->deep_s + tbegin + the_n + d));
// 			// 	printf("\n");

// 			// 	while (prev_reach <= the_n + 1 && the_n >= 0 && the_n+d >= 0 && \
// 			// 		(pattern[the_n] & *(di->deep_s + tbegin + the_n + d)))
// 			// 	{
					
// 			// 		deep_path[d_index] = b;
// 			// 		d_path[d_index++] = mat;
// 		 //            // printf ("mat | ");
// 		 //            // printf ("n: %d d: %d b: %d ", the_n, d, b);
// 		 //            // print_bit4(pattern[the_n]);
// 		 //            // printf ("| ");

// 		 //            --the_n;
// 			// 	}
// 			// }
// 		}

// 		if (d <= 0 && (the_n+1 == -d) && (cur_deep > 0))
// 		{
// 			int cross = 0;
// 			// printf("try cross: from %d \n", b);
// 			int tb = b;
// 			while ((di->all_node[tb]).deep == cur_deep) --tb;

// 			for (;(di->all_node[tb]).deep == cur_deep-1 && cross == 0; --tb)
// 			{
// 				lt = (int)MN(((di->all_node[tb]).len), INT_MAX);
				
// 				for (int td = -lp+1, tdk; td<(int)(MN( e+1, MN(lt, lp))); ++td)
// 				{
// 					tdk = td + lp + 1;
// 					if (the_n + 1 + td == lt && (lv->LV)[tb][ek][tdk] == the_n+1)
// 					{
// 						--cur_deep;
// 						// printf("| cross: from %d to %d\n", b, tb);
// 						d = td;
// 						b = tb;
// 						dk = d + lp + 1;
// 						bn = &(di->all_node[b]);
// 						tbegin = bn->begin;
// 						cross = 1;
// 						if (b > 0 && (di->all_node[b-1]).deep == (di->all_node[b]).deep)
// 							printf("+");
// 						break;
// 					}

// 				}
// 				// printf("%d %d\n", the_n+1, (lv->LV)[tb][ek][dk]);
// 			}
// 			if (cross == 1) continue;
// 		}

// 		--e;
// 		if (e >= 0)
//         {
//             //now ek = e
//             //LV[e][k] is LV[ek][k]
//             --ek;
//             // printf("%d %d %d %d %d | ", e, d, (lv->LV)[b][ek][dk-1], (lv->LV)[b][ek][dk], (lv->LV)[b][ek][dk+1]);
//             if (the_n >= 0 && (lv->LV)[b][ek][dk] == the_n)
//             {
//                 //d = d
// 				deep_path[d_index] = b;
//                 d_path[d_index++] = mis;
//                 --the_n;
//                 // printf ("mis | ");
//             }
//             else if ((lv->LV)[b][ek][dk+1]  >= 0 && (lv->LV)[b][ek][dk+1] >= the_n)
//             {
            	
// 				deep_path[d_index] = b;
//                 d_path[d_index++] = ins;
//                 --the_n;
//                 // printf ("ins | ");
//                 d = d + 1;
//             }
//             else if ((lv->LV)[b][ek][dk-1] >= (the_n - 1))
//             {
// 				deep_path[d_index] = b;
//                 d_path[d_index++] = del;
//                 // printf ("del | ");
//                 d = d - 1;
//             }
//             else
//             {
//                 //reach bound
//                 fprintf(stderr, "retrive cigar error\n");
//                 printf("retrive cigar error\n");
//                 break;
//             }
//         }
// 	}

// 	//short sam format
// 	// puts("");
// 	int o=0;
// 	char tmp = ((d_path[d_index-1] == 'X') ? 'M' : d_path[d_index-1]);
// 	// char tmp = (d_path[d_index-1]);
// 	char c = tmp;
//     for (int i=d_index-1; i>=0; --i)
//     {
//     	tmp = ((d_path[i-1] == 'X') ? 'M' : d_path[i-1]);
//     	++o;
//     	// tmp = (d_path[d_index-1]);
//     	if (tmp!= c)
//     	{
//     		printf("%d%c", o, c);
//     		c = tmp;
//     		o = 0;
//     	}
//         // printf("%c", d_path[i]);
//     }
//     if (o!=0) printf("%d%c", o, c);
//     puts("");

//     // long alignment format
//     if (LOGA)
//     {
// 		puts("");
// 		for (int i=d_index-1; i>=0; --i)
// 	    {
// 	    	printf("%c", d_path[i]);
// 	    }
// 	    puts("");
// 	    for (int i=d_index-1; i>=0; --i)
// 	    {
// 	        printf("%d", deep_path[i]);
// 	    }
// 	    puts("");

// 	    b = -1;
// 	    bn = &(di->all_node[b]);
// 		tbegin = bn->begin;
// 	    for (int i=d_index-1, j; i>=0; --i)
// 	    {
// 	        if (d_path[i] == ins) printf("_");
// 	        else
// 	        {
// 	        	if (deep_path[i] != b)
// 	        	{
// 	        		b = deep_path[i];
// 	        		bn = &(di->all_node[b]);
// 	        		tbegin = bn->begin;
// 	        		j = 0;
// 	        		// printf("[b:%d]", b);
// 	        	}
// 	            int o = *(di->deep_s + tbegin + (j++));
// 	            print_bit4_(o);
// 	        }
// 	    }
// 	    puts("");

// 	    for (int i=d_index-1, j=0; i>=0; --i)
// 	    {
// 	        if (d_path[i] == del) printf("_");
// 	        else 
// 	        {
// 	        	int o = pattern[j++];
// 	            print_bit4_(o);
// 	        }
// 	    }
// 	    puts("");
//     	puts("");
//     }

//     free(d_path);
//     free(deep_path);
// }