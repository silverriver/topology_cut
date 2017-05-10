#pragma once
#include <vector>
#include "geometry.h"
#include <fstream>

namespace TC
{

class block_system
{
public:
	std::vector<point_3> vpo;
	std::vector <plane_3> vpl;
	std::vector <face> vfa;			//原始边界，未砍树之前
	std::vector <block> vblk;		//原始块体，未砍树之前

	std::vector <face> vfa_ext;		//块体的外表面face
	std::vector <block> vblk_ext;	//每个块体的外表面，其中的索引对应数组vfa_ext中的面
	std::vector <block> vblk_int;	//每个块体的内表面

	//记录原始的数据，包括domain和裂隙
	stl_model domain;				
	std::vector <poly_frac> vpf;
	std::vector <disc_frac> vdf;		//圆盘裂隙和多边形裂隙分开处理

public:
	int init_domain (std::ifstream & infile);

	int add_polyf (const poly_frac& pf, const poly_frac::Type &type = poly_frac::Type::Both_side);

	int identify ();

	int disc_frac_parser (std::ifstream &infile, const FT& dx, const int & n = 10);

	int output_ved (std::ofstream &outfile) const;

	int output_fac (std::ofstream &outfile, const int&i) const;

	int output_blk (std::ofstream &outfile, const int &i) const;

	int output_outer_bounds (std::ofstream &outfile, const int &i) const;

	double cal_vol_estimate (const block &b) const;

	FT cal_vol (const block &b) const;
public:
	std::vector <edge> ved;		//每个包含在块体系统内或边界上的edge
	std::vector <std::vector<size_t>> ved_pli;		//每个edge在哪些平面上
	std::vector <std::vector<size_t>> vpl_edi;		//每个平面上有哪些edge
	std::vector <csection> vsec;		//记录每个平面pl与domain的截面
	
	std::vector <csection> vsec_pf;		//记录每个平面pl上有哪些多边形裂隙和domain边界

	//生成各个平面和domian的截面，包括内部的点和边界上的点
	//这些计算出来的截面会在生成edge和生成face的时候用到
	int init_csections (const int& pli, csection &sec);

	//初始化vsec_pf,记得在生成各个face之前调用
	int init_vsec_pf ();

	//计算所有的edge，并且给ved和vpo赋值
	int gen_edges ();

	//Generate the loops on each plane. 
	//pli: index to the plane
	//vl:  results
	int gen_loop (const int& pli, std::vector<loop> &vl);

	//将loop转换为face，并储存起来,vl是split过之后的. pli表示所在平面的索引
	int loop2face (const int& pli, const std::vector<loop> &vl);

	//生成原始块体，也就是未砍树之前的块体，可以应用在渗流程序中
	int gen_blk ();

	//提取块体的外表面边界，也就是给vfa_ext, vblk_ext赋值
	int extract_outer_surfaces();

	//检查数据的协调性
	int check_face_loop_edge () const;

	int gen_outer_boundary ();

	int sort_edge ()
	{
		edge_compare ec;
		std::sort (ved.begin (),ved.end (),ec);
	}

};

}