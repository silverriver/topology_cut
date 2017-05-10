#include "geometry.h"
#include "block_system.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <CGAL\Timer.h>
using namespace std;

namespace TC
{


int block_system::init_domain (std::ifstream & infile) 
{
	domain.clear();
	int res1 = domain.init_stl (infile);
	int res2 = domain.check ();
	if (res1 == 0 && res2 == 0)
	{
		for (int i(0);i!=domain.vt.size ();++i)		//记录每个面的转向
		{
			plane_3 temp (domain.get_plane(i));
			domain.vpli.push_back (add_set_pl(vpl,temp));
			if (temp.orthogonal_vector() * vpl[domain.vpli.back ()].orthogonal_vector() > 0)
				domain.vtri_ori.push_back (true);
			else
				domain.vtri_ori .push_back (false);			
		}
		return 0;
	}
	else
		return 1;
}

int block_system::add_polyf (const poly_frac& pf, const poly_frac::Type &type)
{
	vector_3 normal (oriented_normal (pf.vl[0]));
	if (normal == vector_3 (0,0,0))
		return 1;
	
	vpf.push_back (pf);
	vpf.back ().plane_id = add_set_pl(vpl,plane_3(pf.vl[0][0],normal));
	vpf.back ().frac_id = (int) vpf.size ();
	vpf.back ().type = type;
	vpf.back ().proj_plane(vpl[vpf.back ().plane_id]);		//以防万一

	if (normal * vpl[vpf.back ().plane_id].orthogonal_vector() < 0)
		vpf.back ().revise_loops();		//需要翻转各个loop

	return 0;
}

int block_system::identify ()
{
	ofstream outfile;						//DEBUG

	cout<<"Init. cross section"<<endl;
	vsec.resize (vpl.size ());
	for (int i(0);i!=vpl.size ();++i)
		init_csections(i,vsec[i]);
	cout<<"Init. cross section complete"<<endl;

	cout<<"Init. cross sections for fractures"<<endl;
	init_vsec_pf();
	cout<<"Init. cross sections for fractures complete"<<endl;

	cout<<"Generating edges"<<endl;
	gen_edges();
	cout<<"Edge generation complete"<<endl;

	CGAL::Timer timer;
	//生成各个face
	timer.reset ();
	timer.start();
	for (int i(0);i!=vpl.size ();++i)
	{
		cout<<i<<endl;
		vector<loop> vl;
		gen_loop (i,vl);
		loop2face (i,vl);
	}
	timer.stop();
	cout<<"tracing face:"<<timer.time()<<endl;
	timer.reset ();

	//DEBUG
	outfile.open("face_debug.stl");			//DEBUG
	for (int i(0);i!=vfa.size ();++i)		//DEBUG
		output_fac(outfile,i);				//DEBUG
	outfile.close ();						//DEBUG
	outfile.open ("edge_debug.txt");		//DEBUG
	output_ved(outfile);					//DEBUG
	outfile.close ();						//DEBUG

	timer.start ();
	gen_blk();
	cout<<"tracing blks:"<<timer.time ()<<endl;
	timer.stop (); timer.reset ();
	timer.start ();
	extract_outer_surfaces();
	cout<<"extracting blks' outer surfaces:"<<timer.time ()<<endl;
	return 0;
}

int block_system::disc_frac_parser (std::ifstream &infile, const FT& dx, const int&n)
{
	if (!infile)
		return 1;
	int count(0);
	infile>>count;
	for (int i(0);i!=count;++i)
	{
		cout<<i<<"/"<<count<<endl;
		if (!infile)
			return 0;
		poly_frac temp;
		FT x,y,z,dip_dir,dip,r;
		infile>>x>>y>>z>>dip_dir>>dip>>r>>temp.aperture >>temp.cohesion>>temp.f_angle;
		disc2pf(point_3 (x,y,z),r,plane_3 (point_3(x,y,z),dip2normal(dip_dir,dip,dx)),temp,n);
		add_polyf(temp);
	}
	return 0;
}

int block_system::output_ved (std::ofstream &outfile) const
{
	if (!outfile)
		return 1;
	for (int i(0);i!=ved.size ();++i)
		out_line_3_rhino(outfile,segment_3(vpo[ved[i].n[0]],vpo[ved[i].n[1]]));
	return 0;
}

int block_system::output_fac (std::ofstream &outfile, const int&i) const
{
	stringstream ss;
	ss<<i;
	out_face_STL(outfile,vpo,vfa[i],vpl[vfa[i].pli],true,ss.str());
	return 0;
}

int block_system::output_blk (std::ofstream &outfile, const int &bi) const
{
	stringstream ss;
	ss<<bi;
	outfile<<"solid BLK"<<ss.str()<<endl;
	for (int i(0);i!=vblk[bi].vfi.size ();++i)
	{
		plane_3 pl;
		if (vblk[bi].vfi[i].second )
			pl = vpl[vfa[vblk[bi].vfi[i].first ].pli];
		else
			pl = vpl[vfa[vblk[bi].vfi[i].first ].pli].opposite();

		out_face_STL(outfile,vpo,vfa[vblk[bi].vfi[i].first],pl,false);
	}
	outfile<<"endsolid BLK"<<ss.str ()<<endl;
	return 0;

}

int block_system::output_outer_bounds (std::ofstream &outfile, const int &bi) const
{
	stringstream ss;
	ss<<bi;
	outfile<<"solid BLK"<<ss.str()<<endl;
	for (int i(0);i!=vblk_ext[bi].vfi.size ();++i)
	{
		plane_3 pl;
		if (vblk_ext[bi].vfi[i].second )
			pl = vpl[vfa_ext[vblk_ext[bi].vfi[i].first ].pli];
		else
			pl = vpl[vfa_ext[vblk_ext[bi].vfi[i].first ].pli].opposite();

		out_face_STL(outfile,vpo,vfa_ext[vblk_ext[bi].vfi[i].first],pl,false);
	}
	outfile<<"endsolid BLK"<<ss.str ()<<endl;
	return 0;
}

double block_system::cal_vol_estimate (const block &b) const
{
	using CGAL::to_double;
	double res(0);
	for (int i(0);i!=b.vfi.size ();++i)
	{
		const point_3 &p1 (vpo[vfa[b.vfi[i].first ].vl[0].vpi[0]]);
		for (int j(0);j!=vfa[b.vfi[i].first].vl.size ();++j)
		{
			for (int k(0);k!=vfa[b.vfi[i].first].vl[j].vei.size ();++k)
			{
				// + - -
				// + + +
				// - + -
				// - - +
				edge e (ved[vfa[b.vfi[i].first].vl[j].vei[k]]);
				if (b.vfi[i].second != vfa[b.vfi[i].first].vl[j].ve_ori[k])
					e = e.opposite();
				const point_3 &p2 (vpo[e.n[0]]), &p3 (vpo[e.n[1]]);
				res += to_double (p1.x() * p2.y() * p3.z() + p1.y() * p2.z() * p3.x() + p1.z() * p2.x() * p3.y()
					-p1.z() * p2.y() * p3.x() - p1.y() * p2.x() * p3.z() - p1.x() * p2.z() * p3.y());
			}
		}
	}
	return res/6.0;
}

FT block_system::cal_vol (const block &b) const
{
	FT res(0);
	for (int i(0);i!=b.vfi.size ();++i)
	{
		const point_3 &p1 (vpo[vfa[b.vfi[i].first ].vl[0].vpi[0]]);
		for (int j(0);j!=vfa[b.vfi[i].first].vl.size ();++j)
		{
			for (int k(0);k!=vfa[b.vfi[i].first].vl[j].vei.size ();++k)
			{
				// + - -
				// + + +
				// - + -
				// - - +
				edge e (ved[vfa[b.vfi[i].first].vl[j].vei[k]]);
				if (b.vfi[i].second != vfa[b.vfi[i].first].vl[j].ve_ori[k])
					e = e.opposite();
				const point_3 &p2 (vpo[e.n[0]]), &p3 (vpo[e.n[1]]);
				res += (p1.x() * p2.y() * p3.z() + p1.y() * p2.z() * p3.x() + p1.z() * p2.x() * p3.y()
					-p1.z() * p2.y() * p3.x() - p1.y() * p2.x() * p3.z() - p1.x() * p2.z() * p3.y());
			}
		}
	}
	return res/6.0;
}

int block_system::gen_loop (const int& i, std::vector <loop> &vl)
{
	//vector<loop> vl_edge;	//对于每个loop，记录其边在ved上的索引。vl_edge[i].vpi[i]记录的是边(vl[i].vpi[i], vl[i].vpi[i+1])		【暂时没有用到】
	//vector<vector<bool>> vl_edge_ori;		//记录每个loop的边方向，true：在这个loop中的这个边和vl_edge中索引对应的边方向相同，false：方向相反	【暂时没有用到】

	vector<int> vet (vpl_edi[i].size (),0);		//每个edge有哪些已经被利用了0:都还没有被添加，1：正向已添加，2：反向已添加,3:都已经被添加
	int count(0);			//已经有多少边被添加了
	bool tracing (false);		//是否正在完成一个loop的过程中
	int vpl_edi_size ((int)vpl_edi[i].size ());

	while (count < 2*vpl_edi_size)
	{
		if (!tracing)			//要开始一个新的loop
		{
			vl.push_back (loop());

			for (int j(0);j!=vpl_edi_size;++j)
			{
				bool temp(false);
				switch (vet[j])		//找一个边继续开始
				{
				case 0:				//这个边没有被标记过
					temp = true;
					vl.back().vpi.push_back (ved[vpl_edi[i][j]].n[0]);
					vl.back ().vpi.push_back (ved[vpl_edi[i][j]].n[1]);
					vl.back ().vei.push_back ((int)vpl_edi[i][j]);
					vl.back ().ve_ori.push_back (true);
					vet[j] = 1; count++; tracing= true;
					break;
				case 1:				//这个边正向被标记过
					temp = true;
					vl.back ().vpi.push_back (ved[vpl_edi[i][j]].n[1]);
					vl.back ().vpi .push_back (ved[vpl_edi[i][j]].n[0]);
					vl.back ().vei.push_back ((int)vpl_edi[i][j]);
					vl.back ().ve_ori.push_back (false);
					vet[j] = 3; count++; tracing = true;
					break;
				case 2:				//这个边反向被标记过
					temp = true;
					vl.back().vpi.push_back (ved[vpl_edi[i][j]].n[0]);
					vl.back ().vpi.push_back (ved[vpl_edi[i][j]].n[1]);
					vl.back ().vei.push_back ((int)vpl_edi[i][j]);
					vl.back ().ve_ori.push_back (true);
					vet[j] = 3; count++; tracing= true;
					break;
				default:
					break;
				}
				if (temp)
					break;
			}
		}
		else
		{								//正在完成一个loop的过程中,这个未完成的loop是vl.back()
			vector<int> ve_cand;		//备选边的集合,记录的是vpl_edi[i]中的元素的索引
			vector<bool> ve_dir;
			ve_cand.reserve (10);
			ve_dir.reserve (10);
			for (int j(0);j!=vpl_edi_size;++j)
			{
				if (vet[j] == 3) continue;		//这个边已经被用过了
				if ((vl.back ().vpi.back() == ved[vpl_edi[i][j]].n[0]) && (vet[j] == 0 || vet[j] == 2))
				{
					ve_cand.push_back (j); ve_dir.push_back (true);		//正向还没有添加的时候
				}

				if ((vl.back ().vpi.back() == ved[vpl_edi[i][j]].n[1]) && (vet[j] == 0 || vet[j] == 1))
				{
					ve_cand.push_back (j); ve_dir.push_back (false);
				}
			}
			if (ve_cand.empty())
				throw logic_error ("identify(), face tracing, no candiate edge");
			int maxi (0);		//记录最右的边在ve_cand中的索引
			edge e_prev (edge(vl.back ().vpi[vl.back ().vpi.size ()-2], vl.back().vpi.back ()));
			edge e_next (ved[vpl_edi[i][ve_cand[maxi]]]);		//最右的边
			if (!ve_dir[maxi]) e_next = e_next.opposite();

			right_most_compare_3 compare (vpl[i],vpo,e_prev);
			for (int j(0);j!=ve_cand.size ();++j)
			{
				if (ve_dir[j])		//正在处理一个正向边
				{
					if (compare (e_next,ved[vpl_edi[i][ve_cand[j]]]))
					{
						maxi = j; e_next = ved[vpl_edi[i][ve_cand[j]]];
					}
				}
				else			//正在处理一个反向边
					if (compare (e_next, ved[vpl_edi[i][ve_cand[j]]].opposite()))
					{
						maxi = j; e_next = ved[vpl_edi[i][ve_cand[j]]].opposite();
					}
			}

			switch (vet[ve_cand[maxi]])
			{
			case 0:		//找到的最右边没有被添加过
				if (ve_dir[maxi])
					vet[ve_cand[maxi]] = 1;
				else
					vet[ve_cand[maxi]] = 2;

				break;
			case 1:		//找到的最右边正向被添加过
				if (ve_dir[maxi])
					throw logic_error ("identify(), face tracing edge marked again");
				else
					vet[ve_cand[maxi]] = 3;
				break;
			case 2:		//找到的最右边反向被添加过
				if (ve_dir[maxi])
					vet[ve_cand[maxi]] = 3;
				else
					throw logic_error ("identify(), face tracing edge marked again");
				break;
			case 3:
				throw logic_error ("identify(), face tracing edge marked again");
			default:	break;
			}
			count ++;
			if (vl.back ().vpi [0] == e_next.n[1])		//如果要添加的边等于这个loop的第一个边，则跳出
				tracing = false;
			else
				vl.back ().vpi.push_back (e_next.n[1]);

			vl.back ().vei.push_back ((int)vpl_edi[i][ve_cand[maxi]]);
			vl.back ().ve_ori.push_back (ve_dir[maxi]);
		}
		
	}
	if (tracing)
		throw logic_error ("trac_face(), edge is consumed before finishing the loop");

	vector<loop> vl_temp;
	vl_temp.reserve (vl.size ()*5);
	for (int j(0);j!= vl.size ();++j)
		split_loop(vl[j],vl_temp);
	vl.swap (vl_temp);

	return 0;
}

int block_system::loop2face (const int& pli,const std::vector<loop> &vl)
{
	vector<int> vlt (vl.size ());		//每个loop的type,0:退化边界，1：顺时针旋转，2：逆时针旋转
	vector<vector<int>> vl_containing (vl.size ());		//每个loop包含其他loop的索引

	for (int i(0);i!=vl.size ();++i)
	{
		switch (vl[i].vpi.size ())
		{
		case 0:
		case 1:
			throw logic_error("loop2face, illegal loops"); break;
		case 2:
			vlt[i] = 0; break;
		default:
			if (oriented_normal(vpo,vl[i].vpi) * vpl[pli].orthogonal_vector() > 0)
				vlt[i] = 2;
			else
				vlt[i] = 1;
			break;
		}
	}

	vector<index_l> temp_vl (vl.size ());		//将vl中顺时针旋转的loop都改成逆时针的,并且只复制vpi部分
	for (int i(0);i!=vl.size();++i)
	{
		if (vlt[i] == 0)
			continue;
		else if (vlt[i] == 1)
		{
			temp_vl[i] = vl[i].vpi;
			reverse (temp_vl[i].begin (), temp_vl[i].end ());
		}
		else 
			temp_vl[i] = vl[i].vpi;
	}

	//判断loop之间的包含关系
	for (int i(0);i!=vl.size ();++i)
	{
		if (vlt[i] == 0)
			continue;

		for (int j(0);j!=vl.size ();++j)
		{
			if (i==j) continue;
			if (vlt[j] != 0)
			{
				//如果这两个loop相等就不认为是包含关系，何必互相伤害
				if ((temp_vl[i] != temp_vl[j]) && (polygon_in_polygon_3(vpo,temp_vl[i],temp_vl[j],vpl[pli]) == 1))
					vl_containing[i].push_back (j);
			}
			else
			{
				switch (point_in_polygon_3(vpo,temp_vl[i],vpl[pli],CGAL::midpoint(vpo[vl[j].vpi[0]], vpo[vl[j].vpi[1]])))
				{
				case 0:
					throw logic_error ("loop2face(), edge on loop"); break;
				case 1:
					vl_containing[i].push_back (j);
				} 
			}
		}
	}

	vector<int> vl_fai (vl.size (),-1);		//记录vl中每个outer loop被添加到哪个face中了
	for (int i(0);i!=vl.size ();++i)
		if (vlt[i] == 2)
		{
			for (int j(0);j<vl_containing[i].size ();++j)
			{
				int temp_li (vl_containing[i][j]);
				if (vlt[temp_li] == 2)
				{
					remove_set (vl_containing[i],vl_containing[temp_li]);
					remove_set (vl_containing[i],temp_li);
					--j;
				}
			}
		}

	vector<face> temp_vfa;
	temp_vfa.reserve (vl.size ());
	for (int i(0);i!=vl.size ();++i)
	{
		if (vlt[i] == 2)
		{
			temp_vfa.push_back(face());
			temp_vfa.back ().vl.push_back (vl[i]);
			temp_vfa.back ().pli = pli;
			if (vl_fai[i] == -1)
				vl_fai[i] = (int)temp_vfa.size ()-1;
			else
				throw logic_error ("loop2face, loop has been added");
			for (int j(0);j!=vl_containing[i].size ();++j)
			{
				temp_vfa.back ().vl.push_back (vl[vl_containing[i][j]]);
				if (vl_fai[vl_containing[i][j]] == -1)
					vl_fai[vl_containing[i][j]] = (int)temp_vfa.size ()-1;
				else
					throw logic_error ("loop2face, loop has been added");
			}
		}
	}

	for (int i(0);i!=temp_vfa.size ();++i)
	{
		point_3 temp (get_point_in_face_3(vpo,vpl[pli],temp_vfa[i]));
		//如果这个face上的点即在横截面上，也在裂隙面上
		if ((point_in_polygon_3(vsec[pli].vp,vsec[pli].vl,vpl[pli],temp) != 2)
			&& (point_in_polygon_3(vsec_pf[pli].vp,vsec_pf[pli].vl,vpl[pli],temp) !=2))
			vfa.push_back (temp_vfa[i]);
	}
	return 0;
}

int block_system::gen_edges ()
{
	std::vector <segment_3> vseg;
	std::vector <point_3> vtempp;
	std::vector <bbox_3> vbox_pf;		//bbox for each poly frac
	std::vector <bbox_3> vbox_seg;		//bbox for each segment in vseg
	std::vector<int> vseg_type;			//whether the segment in vseg is in the domain
										//0: Not judged yet. 1:inside. 2:outside

	CGAL::Timer timer;
	cout<<"Initializing bbox ";
	timer.start ();
	//Initialize bbox for each poly fracture
	vbox_pf.reserve (vpf.size ());
	for (int i(0);i!=vpf.size ();++i)
		vbox_pf.push_back (vpf[i].init_bbox());

	vseg.reserve (domain.vt.size ()*3 + vpf.size ()*30);
	vbox_seg.reserve (domain.vt.size ()*3 + vpf.size ()*30);
	for (int i(0);i!=domain.vt.size();++i)
	{
		add_set_seg_3(vseg,segment_3 (domain.vp[domain.vt[i].n[0]],domain.vp[domain.vt[i].n[1]]),vbox_seg);
		add_set_seg_3(vseg,segment_3 (domain.vp[domain.vt[i].n[0]],domain.vp[domain.vt[i].n[2]]),vbox_seg);
		add_set_seg_3(vseg,segment_3 (domain.vp[domain.vt[i].n[2]],domain.vp[domain.vt[i].n[1]]),vbox_seg);
	}

	for (int i(0);i!=vpf.size ();++i)
		for (int j(0);j!=vpf[i].vl.size ();++j)
			for (int k(0);k!=vpf[i].vl[j].size ();++k)
			{
				int k1 (k+1);
				if (k1 >= vpf[i].vl[j].size ())
					k1 = 0;
				add_set_seg_3 (vseg, segment_3 (vpf[i].vl[j][k], vpf[i].vl[j][k1]), vbox_seg);
			}

	timer.stop ();
	cout<<timer.time ()<<endl;

	timer.reset ();
	timer.start ();
	cout<<"Generating all the intersected segments "; 
	//计算所有的交线段
	for (int i(0);i!= vpf.size ();++i)
	{
		cout<<i<<endl;
		for (int j (i+1); j<vpf.size ();++j)
		{
			if (!CGAL::do_overlap (vbox_pf[i],vbox_pf[j]))
				continue;
			intersection(vpf[i], vpl[vpf[i].plane_id ],
				vpf[j], vpl[vpf[j].plane_id],
				vseg,vtempp,vbox_seg);
		}
	}
	
	for (int i(0);i!= domain.vt.size ();++i)
	{
		cout<<i<<"/"<<domain.vt.size ()<<endl;
		plane_3 temp_pl(vpl[domain.vpli[i]]);
		if (!domain.vtri_ori[i])
			temp_pl = temp_pl.opposite();
		poly_frac pf;

		bbox_3 df_bbox (domain.init_bbox (i));
		
		pf.vl.push_back (poly_frac::f_loop());
		pf.vl.back ().push_back (domain.vp[domain.vt[i].n[0]]);
		pf.vl.back ().push_back (domain.vp[domain.vt[i].n[1]]);
		pf.vl.back ().push_back (domain.vp[domain.vt[i].n[2]]);
		for (int j(0);j!=vpf.size ();++j)
		{
			if (!CGAL::do_overlap(df_bbox,vbox_pf[j]))
				continue;
			intersection (pf,temp_pl,vpf[j],vpl[vpf[j].plane_id],
						vseg,vtempp,vbox_seg);
		}
	}
	
	//用vtempp中的点分裂各个线段
	sort(vtempp.begin (),vtempp.end (), point_compare_3());
	vtempp.resize (unique (vtempp.begin (),vtempp.end ()) - vtempp.begin ());
	for (int i(0);i!=vtempp.size ();++i)
		split_seg(vseg, vbox_seg, vtempp[i]);

	timer.stop ();
	cout<<timer.time ()<<endl;
	timer.reset ();
	
	

	//剔除在研究区域以外的线段
	timer.start ();
	cout<<"Deleting segments out of domain ";
	vseg_type.resize (vseg.size (),0);					//每个线段是否应该进入下一步的分析,0:还没有判断，1：要进入，2：不要进入
	vector<int> vpl_seg_count (vpl.size (),0);			//每个平面上有多少线段已经被判断过了
	vector<vector<size_t>> vseg_pli (vseg.size ());		//每个线段在哪些平面上
	vector<vector<size_t>> vpl_segi (vpl.size ());		//每个平面上的线段
	
	//vector<vector<size_t>> vpl_domainti (vpl.size ());	//每个平面上有多少domain的三角面
	//for (int i(0);i!=domain.vpli.size ();++i)
	//{
	//	vpl_domainti[domain.vpli[i]].push_back (i);
	//}

	for (int i(0);i!=vseg.size ();++i)
	{
		for (int j(0);j!=vpl.size ();++j)
		{
			if (vpl[j].has_on(vseg[i].target()) && vpl[j].has_on(vseg[i].source ()))
			{
				vseg_pli[i].push_back (j);
				vpl_segi[j].push_back (i);
			}
		}
	}

	for (int i(0);i!=vpl.size ();++i)
	{
		cout<<i<<endl;

		if (vpl_seg_count[i] == vpl_segi[i].size ())
			continue;		//如果每一个线段都被判断了直接继续

		for (int j(0);j!=vpl_segi[i].size ();++j)
		{
			if (vseg_type[vpl_segi[i][j]] != 0)
				continue;		//这个线段已经被判断过了
			vseg_type[vpl_segi[i][j]] = 1;	//默认要进入
			vpl_seg_count[i] ++;
			
			point_3 mid(CGAL::midpoint(vseg[vpl_segi[i][j]].source (), vseg[vpl_segi[i][j]].target ()));
			if (point_in_polygon_3(vsec[i].vp,vsec[i].vl,vpl[i],mid) == 2)
				vseg_type [vpl_segi[i][j]] = 2;
		}
	}
	timer.stop ();
	cout<<timer.time()<<endl;
	timer.reset ();

	//将包含在研究区域内的边都添加到ved中
	timer.start ();
	cout<<"Generating ved ";
	vpl_edi.resize (vpl.size ());
	ved_pli.reserve (vseg.size ());
	for (int i(0);i!=vseg.size ();++i)
	{
		switch (vseg_type[i])
		{
		case 0:
			throw logic_error ("gen_edges, edge unmarked"); break;
		case 1:
			ved.push_back (edge(add_set (vpo,vseg[i].source ()), add_set(vpo,vseg[i].target ())));
			ved.back ().sort ();
			ved_pli.push_back (vseg_pli[i]);
			for (int j(0);j!=vseg_pli[i].size ();++j)
				vpl_edi[vseg_pli[i][j]].push_back ((int)ved.size ()-1);
			break;
		default:
			break;
		}
	}
	timer.stop ();
	cout<<timer.time ()<<endl;
	return 0;
}

int block_system::init_csections (const int& pli, csection &sec)
{
	//可以优化一下的，如果直接计算所有的截面的话只需要遍历一遍domain.vt就可以了
	sec.clear ();
	domain.plane_intersect(vpl[pli],sec);
	for (int i(0);i!=domain.vt.size ();++i)
	{
		if (domain.vpli[i] != pli)
			continue;
		sec.vl.push_back (index_l());
		sec.vp.push_back (domain.vp[domain.vt[i].n[0]]);
		sec.vp.push_back (domain.vp[domain.vt[i].n[1]]);
		sec.vp.push_back (domain.vp[domain.vt[i].n[2]]);
		if (domain.vtri_ori[i])
		{
			sec.vl.back ().push_back ((int)sec.vp.size () - 3);
			sec.vl.back ().push_back ((int)sec.vp.size () - 2);
			sec.vl.back ().push_back ((int)sec.vp.size () - 1);
		}
		else
		{
			sec.vl.back ().push_back ((int)sec.vp.size () - 1);
			sec.vl.back ().push_back ((int)sec.vp.size () - 2);
			sec.vl.back ().push_back ((int)sec.vp.size () - 3);
		}
	}
	return 0;
}

int block_system::init_vsec_pf ()
{
	std::vector <std::vector <size_t>> vpl_dfai;		//记录每个平面上有哪些domain的三角面索引
	std::vector <std::vector <size_t>> vpl_pfi;			//记录每个平面上有哪些多边形裂隙
	
	vsec_pf.resize (vpl.size ());
	vpl_dfai.resize (vpl.size ());
	vpl_pfi.resize (vpl.size ());

	for (int i(0);i!=domain.vt.size ();++i)
		vpl_dfai[domain.vpli[i]].push_back (i);
	for (int i(0);i!=vpf.size ();++i)
		vpl_pfi[vpf[i].plane_id].push_back (i);

	for (int i(0);i!=vpl.size ();++i)
	{
		vsec_pf[i].clear ();
		for (int j(0);j!=vpl_dfai[i].size ();++j)
		{
			vsec_pf[i].vl.push_back (index_l());
			vsec_pf[i].vp.push_back (domain.vp[domain.vt[vpl_dfai[i][j]].n[0]]);
			vsec_pf[i].vp.push_back (domain.vp[domain.vt[vpl_dfai[i][j]].n[1]]);
			vsec_pf[i].vp.push_back (domain.vp[domain.vt[vpl_dfai[i][j]].n[2]]);
			if (domain.vtri_ori[vpl_dfai[i][j]])
			{
				vsec_pf[i].vl.back ().push_back ((int) vsec_pf[i].vp.size ()-3);
				vsec_pf[i].vl.back ().push_back ((int) vsec_pf[i].vp.size ()-2);
				vsec_pf[i].vl.back ().push_back ((int) vsec_pf[i].vp.size ()-1);
			}
			else
			{
				vsec_pf[i].vl.back ().push_back ((int) vsec_pf[i].vp.size ()-1);
				vsec_pf[i].vl.back ().push_back ((int) vsec_pf[i].vp.size ()-2);
				vsec_pf[i].vl.back ().push_back ((int) vsec_pf[i].vp.size ()-3);
			}
		}
		for (int j(0);j!=vpl_pfi[i].size ();++j)
			for (int k(0);k!= vpf[vpl_pfi[i][j]].vl.size ();++k)
			{
				vsec_pf[i].vl.push_back (index_l());
				for (int m(0);m!=vpf[vpl_pfi[i][j]].vl[k].size ();++m)
				{
					vsec_pf[i].vl.back ().push_back ((int) vsec_pf[i].vp.size ());
					vsec_pf[i].vp.push_back (vpf[vpl_pfi[i][j]].vl[k][m]);
				}
			}
	}
	return 0;
}

int block_system::check_face_loop_edge () const
{
	std::vector <std::vector<size_t>> temp_ved_pli(ved.size ());		//每个edge在哪些平面上
	std::vector <std::vector<size_t>> temp_vpl_edi(vpl.size ());		//每个平面上有哪些edge
	//检查ved
	for (int i(0);i!=ved.size ();++i)
	{
		if (ved[i].n[0] >= ved[i].n[1])
			return 1;
		for (int j(0);j!=vpl.size ();++j)
			if (vpl[j].has_on (vpo[ved[i].n[0]]) && vpl[j].has_on (vpo[ved[i].n[1]]))
			{
				temp_ved_pli[i].push_back (j);
				temp_vpl_edi[j].push_back (i);
			}
	}

	//检查ved_pli
	for (int i(0);i!=ved.size ();++i)
	{
		//ved_pli[i] is included in temp_ved_pli[i]
		for (int j(0);j!=ved_pli[i].size ();++j)
			if (find (temp_ved_pli[i].begin (), temp_ved_pli[i].end (),ved_pli[i][j]) == temp_ved_pli[i].end ())
				return 2;

		for (int j(0);j!=temp_ved_pli[i].size ();++j)
			if (find (ved_pli[i].begin (), ved_pli[i].end (),temp_ved_pli[i][j]) == ved_pli[i].end ())
				return 3;
	}

	//检查vpl_edi
	for (int i(0);i!=vpl.size ();++i)
	{

		for (int j(0);j!=vpl_edi[i].size ();++j)
			if (find (temp_vpl_edi[i].begin (), temp_vpl_edi[i].end (),vpl_edi[i][j]) == temp_vpl_edi[i].end ())
				return 4;

		for (int j(0);j!=temp_vpl_edi[i].size ();++j)
			if (find (vpl_edi[i].begin (), vpl_edi[i].end (),temp_vpl_edi[i][j]) == vpl_edi[i].end ())
				return 5;
	}

	//检查各个face
	for (int i(0);i!=vfa.size ();++i)
	{
		for (int j(0);j!=vfa[i].vl.size ();++j)
		{
			const loop & l = vfa[i].vl[j];
			if (l.vpi.size () != l.vei.size ())
				return 6;
			if (l.vpi.size () != l.ve_ori.size ())
				return 7;
			for (int k(0);k!=l.vpi.size ();++k)
			{
				int k1 (k+1);
				if (k1 >= l.vpi.size ()) k1=0;
				if (ved[l.vei[k]] != edge(l.vpi[k], l.vpi[k1]))
					return 8;
				if (l.ve_ori[k] && (ved[l.vei[k]].n[0] != l.vpi[k]))
					return 9;
				if ((!l.ve_ori[k]) && (ved[l.vei[k]].n[0] != l.vpi[k1]))
					return 10;
				if (find (ved_pli[l.vei[k]].begin (),ved_pli[l.vei[k]].end (), vfa[i].pli) == ved_pli[l.vei[k]].end ())
					return 11;		//该面的索引没有包含在这个面的edge中
			}
		}
	}
	return 0;
}

int block_system::gen_blk ()
{
	vector<vector<pair<int,bool>>> ved_fai(ved.size ());			//每个edge所在的平面索引和该edge在这个平面的朝向

	for (int i(0);i!=vfa.size ();++i)
		for (int j(0);j!=vfa[i].vl.size ();++j)
		{
			const loop & l = vfa[i].vl[j];
			for (int k(0);k!=l.vei.size ();++k)
				ved_fai[l.vei[k]].push_back (pair<int,bool>(i,l.ve_ori[k]));
		}

	vector<int> vft (vfa.size ());		//每个face有哪些已经被利用了，0：都还没有被添加，1：正向已添加，2：反向已添加,3:都已经被添加
	int count(0);						//有多少face已经被添加了
	bool tracing (false);
	int vfa_size ((int) vfa.size ());

	while (count < 2 * vfa_size)
	{
		if (!tracing)			//要开始一个新的块体
		{
			vblk.push_back (block());

			for (int i(0);i!=vfa_size;++i)
			{
				bool temp (false);
				switch (vft[i])
				{
				case 0:
					temp = true;
					vblk.back ().vfi.push_back (pair<int,bool> (i,true));
					tracing = true;
					break;
				case 1:
					temp = true;
					vblk.back ().vfi.push_back (pair<int,bool> (i,false));
					tracing = true;
					break;
				case 2:
					temp = true;
					vblk.back ().vfi.push_back (pair<int,bool> (i,true));
					tracing = true;
					break;
				default:
					break;
				}
				if (temp)
					break;
			}
		}
		else					//正在完成一个块体
		{
			for (int i(0);i!=vblk.back ().vfi.size ();++i)
			{
				vector_3 fa_n (vpl[vfa[vblk.back ().vfi[i].first].pli].orthogonal_vector());	//当前face的法向向量
				if (!vblk.back ().vfi[i].second)
					fa_n = -fa_n;
				for (int j(0);j!=vfa[vblk.back ().vfi[i].first].vl.size ();++j)
				{
					const loop & l(vfa[vblk.back().vfi[i].first].vl[j]);
					for (int k(0);k!=l.vei.size ();++k)
					{
						//对这个face中包含的每个edge
						vector<pair<int,bool>> cand_fa;			//记录这个当前edge的twin所在的候选面
						bool true_ori_curr (!(l.ve_ori[k] ^ vblk.back ().vfi[i].second));		//当前edge的真正方向，如果当前面是反过来的，则其所有的edge都应该反过来
						for (int m(0);m!=ved_fai[l.vei[k]].size ();++m)
						{
							//if (ved_fai[l.vei[k]][m].first == vblk.back ().vfi[i].first)	//当前面和所处理面是同一个面
							//{
							//	switch (vft[ved_fai[l.vei[k]][m].first])		//如果当前面的twin还没有被别的块体用到，则添加为候选面
							//	{
							//	case 0:
							//		cand_fa.push_back (pair<int,bool>(ved_fai[l.vei[k]][m].first, !vblk.back ().vfi[i].second ));
							//		break;		//只有正反面都没有被用到的时候才可能被添加
							//	case 1:
							//		if (vblk.back ().vfi[i].second)		//当前面和其所在平面同向，但此时正向已被使用，则不可能
							//			throw logic_error("gen_blk, positive face has been used in other blocks");
							//		break;
							//	case 2:
							//		if (!vblk.back ().vfi[i].second )	//当前面和其所在平面反向，但此时反向已被使用，则不可能
							//			throw logic_error("gen_blk, negative face has been used in other blocks");
							//		break;
							//	default:
							//		throw logic_error ("gen_blk, face has been used in other blocks");		//正反面都被用了，理论上这个面就不该被添加
							//		break;
							//	}
							//	continue;
							//}

							switch (vft[ved_fai[l.vei[k]][m].first])
							{
							case 0:		//待处理面两边都没有被利用
								if (true_ori_curr == ved_fai[l.vei[k]][m].second )		
									cand_fa.push_back (pair<int,bool> (ved_fai[l.vei[k]][m].first ,false));//当前边的twin在待处理面的反面
								else
									cand_fa.push_back (pair<int,bool> (ved_fai[l.vei[k]][m].first ,true));//当前边的twin在待处理面的正面
								break;
							case 1:		//待处理面的正向已经被利用，顶多利用其反面
								if (true_ori_curr == ved_fai[l.vei[k]][m].second)
									cand_fa.push_back (pair<int,bool> (ved_fai[l.vei[k]][m].first ,false));//当前边的twin在待处理面的反面
								break;
							case 2:		//待处理面的反向已经被利用，顶多利用其正面
								if (true_ori_curr != ved_fai[l.vei[k]][m].second)
									cand_fa.push_back (pair<int,bool> (ved_fai[l.vei[k]][m].first ,true));//当前边的twin在待处理面的正面
							default:
								break;
							}
						}
						if (cand_fa.empty())
						{
							//if (l.vpi.size () <=2)		//这个loop是退化的loop，所以这种情况有情可原，可能被其他边界用完了
							//	continue;
							//else
								throw logic_error ("gen_blk(), no candidate face");
						}

						//寻找候选平面
						vector<vector_3> cand_dir;				//记录候选面的外表面法向
						cand_dir.reserve (cand_fa.size ());
						for (int m(0);m!=cand_fa.size ();++m)
						{
							cand_dir.push_back (vpl[vfa[cand_fa[m].first].pli].orthogonal_vector());
							if (!cand_fa[m].second)
								cand_dir.back () = - cand_dir.back ();
						}
						vector_3 ev (vector_3 (vpo[ved[l.vei[k]].n[0]], vpo[ved[l.vei[k]].n[1]]));		//当前边的方向向量
						if (!true_ori_curr)
							ev = -ev;
						right_most_compare_3_vector compare(ev,fa_n);
						int maxi(0);
						for (int m(0);m!=cand_fa.size ();++m)
							if (compare (cand_dir[maxi], cand_dir[m]))
								maxi = m;

						//找到了候选平面,直接添加到块体中，重复的直接忽略
						add_set (vblk.back ().vfi,cand_fa[maxi]);
					}
				}
			}

			//完成了这个块体的搜索
			tracing = false;
			count += (int) vblk.back ().vfi.size ();
			for (int i(0);i!=vblk.back ().vfi.size ();++i)
			{
				switch (vft[vblk.back ().vfi[i].first])
				{
				case 0:
					if (vblk.back ().vfi[i].second )
						vft[vblk.back ().vfi[i].first ] = 1;
					else
						vft[vblk.back ().vfi[i].first ] = 2;
					break;
				case 1:
					if (vblk.back ().vfi[i].second )
						throw logic_error ("gen_blk(), face already used");
					else
						vft[vblk.back ().vfi[i].first ] = 3;
					break;
				case 2:
					if (vblk.back ().vfi[i].second )
						vft[vblk.back ().vfi[i].first ] = 3;
					else
						throw logic_error ("gen_blk(), face already used");
					break;
				default:
					break;
				}
			}

		}
	}
	if (tracing)
		throw logic_error ("gen_blk(), face is consumed befor finishing the blocks");
	
	//将各个块体中的退化面挑出来
	int vblk_size ((int)vblk.size ());
	for (int i(0);i!=vblk_size;++i)
	{
		if (vblk[i].vfi.size () <2)
			throw logic_error ("gen_blk, block contain less than 2 faces");
		if (vblk[i].vfi.size () == 2)
			if (vblk[i].vfi[0].first != vblk[i].vfi[1].first)
				throw logic_error ("gen_blk, illegal blk");
		vector<int> v_flag (vblk[i].vfi.size (),0);		//标记是否配对
		int max (1);
		
		for (int j(0);j!=vblk[i].vfi.size ();++j)
		{
			if (v_flag[j] != 0)
				continue;		//已经配对过了
			int k(j+1);
			for (;k< vblk[i].vfi.size ();++k)
				if (vblk[i].vfi[j].first == vblk[i].vfi[k].first )
					break;
			if (k < vblk[i].vfi.size ())	//找到了配对
			{
				v_flag[j] = v_flag[k] = ++max;
				vblk.push_back (block());
				vblk.back ().vfi.push_back (vblk[i].vfi[j]);
				vblk.back ().vfi.push_back (vblk[i].vfi[k]);
			}
		}
		if (max == 1)
			continue;		//没有退化边界
		block temp_b;
		temp_b.vfi .reserve (vblk[i].vfi .size ());
		for (int j(0);j!=vblk[i].vfi.size ();++j)
		{
			if (v_flag[j] != 0)
				continue;
			temp_b.vfi.push_back (vblk[i].vfi[j]);
		}
		if (temp_b.vfi.size () == 0)
		{
			vblk[i] = vblk.back ();
			vblk.pop_back();
		}
		else if (temp_b.vfi.size () < 4)
			throw logic_error ("gen_blk, blk contain less than 4 faces");
		else
			vblk[i] = temp_b;
	}
	return 0;
}

int block_system::extract_outer_surfaces ()
{
	vfa_ext.clear ();
	vblk_ext.clear ();
	vblk_int.clear ();

	vfa_ext.reserve (vfa.size ());
	vblk_ext.reserve (vblk.size ());
	vblk_int.reserve (vblk.size ());
	
	for (int i(0);i!=vblk.size ();++i)
	{
		block &b (vblk[i]);
		switch (b.vfi.size ())
		{
		case 0:
		case 1:
		case 3:	//这些情况下不可能围成一个块体
			throw logic_error("block_system::extract_outer_surfaces, illegal blks");
			break;
		case 2:
			continue;
		default:
			//这种情况块体b中包含四个以上的面
			vblk_ext.push_back (block());
			for (int j(0);j!=b.vfi.size ();++j)
			{
				vfa_ext.push_back (vfa[b.vfi[j].first]);
				vblk_ext.back ().vfi.push_back (pair<int,bool> ((int)vfa_ext.size ()-1, b.vfi[j].second));
			}
			break;
		}
	}
}

int block_system::gen_outer_boundary ()
{
	//确定各个块体的包含关系
	vector<FT> v_vol (vblk.size ());
	vector<point_3> v_p (vblk.size ());		//记录一个在这个块体内的点
	vector<vector<bbox_3>> vvbox(vblk.size ());
	for (int i(0);i!=vblk.size ();++i)
	{
		v_vol[i] = cal_vol(vblk[i]);
		get_vbox_of_blk(vpo,vfa,vblk[i],vvbox[i]);
	}

	vblk_ext.clear ();
	vblk_ext.reserve (vblk.size ());
	vector<vector<int>> v_blk_contain (vblk.size ());	//记录每个块体包含多少个其他块体
	for (int i(0);i!=vblk.size ();++i)
	{
		if (v_vol[i] <= 0)
			continue;
		vblk_ext.push_back (vblk[i]);		//每个体积大于0的都应该是某个块体的外表面
		for (int j(0); j != vblk.size ();++j)
		{
			if (i == j || abs(v_vol[i]) == abs(v_vol[j]))		//体积相等的块体不考虑谁包含谁
				continue;
			if (find (v_blk_contain[j].begin (), v_blk_contain[j].end (),i) != v_blk_contain[j].end ())
				continue;		//块体i已经包含在块体j中了
			if (v_vol[j] != 0)
			{
				if (polyhedron_in_polyhedron_3(vpo,vfa,vpl,vblk[i],vblk[j],v_vol[i],v_vol[j],vvbox[i],vvbox[j]) == 1)
					v_blk_contain[i].push_back (j);
			}
			else
			{
				//TODO blk_j的体积为0的情况下怎么处理
			}
		}
	}


	//TODO 2016-11-14
	


}

}