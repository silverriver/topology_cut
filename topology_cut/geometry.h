#pragma once
#include <vector>
#include <cmath>
#include <utility>
#include <CGAL\Exact_predicates_exact_constructions_kernel.h>
#include <iostream>
#include <exception>
#include <algorithm>
#include <CGAL\Polygon_2.h>
#include <CGAL\Polygon_with_holes_2.h>
#include <CGAL\Polygon_set_2.h>
#include <CGAL\Bbox_3.h>
#include <CGAL\Bbox_2.h>

#define M_PI_2 6.28318530717958647693
#define M_PI 3.14159265358979323846

namespace TC		//short for topology cut
{

typedef CGAL::Exact_predicates_exact_constructions_kernel	K;
typedef K::Point_2		point_2;
typedef K::Point_3		point_3;
typedef K::Plane_3		plane_3;
typedef K::Vector_2		vector_2;
typedef K::Vector_3		vector_3;
typedef K::Segment_2	segment_2;
typedef K::Segment_3	segment_3;
typedef K::Triangle_3	triangle_3;
typedef K::Line_3		line_3;
typedef K::FT			FT;
typedef CGAL::Bbox_3	bbox_3;
typedef CGAL::Bbox_2	bbox_2;
typedef CGAL::Polygon_2<K>	polygon_2;
typedef CGAL::Polygon_with_holes_2<K>	polygon_with_holes_2;
typedef CGAL::Polygon_set_2<K>		polygon_set_2;

typedef std::vector<int> index_l;
bool operator == (const index_l &l1, const index_l &l2);

inline bool operator != (const index_l &l1, const index_l &l2)
{
	return !(l1==l2);
}

inline bool is_simple (const index_l &l)
{
	for (int i(0);i!=l.size ();++i)
		for (int j(i+1);j<l.size ();++j)
			if (l[i] == l[j])
				return false;
	return true;
}

struct loop
{
	index_l vpi;		//indices of points.vpi[0]不一定等于vpi[n]
	index_l vei;		//vei[i] 表示边(vpi[i], vpi[i+1])的索引
	std::vector <bool> ve_ori;	//为每个边分配一个flag，表示这个边的方向。ve_ori[i]=true: loop中这个边与边vei[i]方向相同
	void reverse () {std::reverse (vpi.begin (),vpi.end ());std::reverse (vei.begin (),vei.end ()); std::reverse (ve_ori.begin (),ve_ori.end ());}
};

struct edge
{
	int n[2];
	bool operator== (const edge&e) const
	{
		if ((n[0] == e.n[0] && n[1] == e.n[1]) ||
			(n[1] == e.n[0] && n[0] == e.n[1]))
			return true;
		else
			return false;
	}
	bool operator != (const edge &e) const {return !(this->operator==(e));}
	edge () {}
	edge (const int&n1, const int &n2) {n[0] = n1; n[1] = n2;}
	edge opposite () const {return edge(n[1],n[0]);}
	void sort () {if (n[0] > n[1]) {int temp(n[1]); n[1] = n[0]; n[0] = temp;}}
};

struct face
{
	std::vector<loop> vl;		//indices of loops, vli[0]是外圈，其余为内圈，从pli的正向往下看，vli[0]应该是逆时针旋转，其余loop应该是顺时针旋转

	//int vbi[2];		//index of the block that this face belongs to. a face belongs to at most 2 blocks.
	//bool vft [2];	//vft[i] = true:  the outer normal of this face on block vbi[i] is the positive direction of the plane that this face belongs to.
					//vft[i] = false: the outer normal of this face on block vbi[i] is the negative direction of the plane that this face belongs to.
					//vft[i] = false: vli中各个loop的实际顺序应该翻转一下
	
	int pli;		//the plane holding this face
	face () { /*vbi[0]=vbi[1]=-1; */pli=-1;}
	void revise_loops ()
	{
		for (int i(0);i!=vl.size ();++i)
			std::reverse(vl[i].vpi.begin (),vl[i].vpi.end ());
	}
	bbox_3 init_bbox (const std::vector<point_3> &vp, double tol = 0.000001) const
	{
		using CGAL::to_double;
		FT minx ((vp[vl[0].vpi[0]].x())),maxx((vp[vl[0].vpi[0]].x()));
		FT miny ((vp[vl[0].vpi[0]].y())),maxy((vp[vl[0].vpi[0]].y()));
		FT minz ((vp[vl[0].vpi[0]].z())),maxz((vp[vl[0].vpi[0]].z()));
		for (int i(0);i!=vl.size ();++i)
		{
			for (int j(0);j!=vl[i].vpi.size ();++j)
			{
				const point_3 &p (vp[vl[i].vpi[j]]);
				if (p.x () > maxx) maxx = p.x();
				if (p.x () < minx) minx = p.x();
				if (p.y() >maxy) maxy = p.y();
				if (p.y() <miny) miny = p.y();
				if (p.z() >maxz) maxz = p.z();
				if (p.z() <minz) minz = p.z();
			}
		}
		return bbox_3 (to_double(minx)-tol, to_double(miny) - tol, to_double(minz) -tol,
			to_double (maxx)+tol, to_double (maxy) + tol, to_double (maxz) + tol);
	}
	int init_face (const std::vector<edge> &ve);		//给vl_edge和vl_edge_ori赋值，trace block 之前做一些初始化
	//void merge_loops ();
};

struct block
{
	std::vector <std::pair <int,bool>> vfi;		//indices of faces and direction of each face (true: the outer normal of this face on this block is the positive direction of plane that this face belongs to)
};


//记录一个截面
struct csection
{
	std::vector<point_3> vp;
	std::vector<index_l> vl;
	void clear () {vp.clear(); vl.clear();}
};

struct stl_model
{
	struct tri_3
	{
		int n[3];
		tri_3 () {}
		tri_3 (const int& n1, const int& n2, const int&n3) {n[0]=n1; n[1]=n2; n[2]=n3;}
		int & operator[] (const int&i) 
		{
			if (i<0 || i>=3)
				throw std::range_error ("stl_model::tri_3()");
			return n[i];
		}
		const int & operator[] (const int&i) const 
		{
			if (i<0 || i>=3)
				throw std::range_error ("stl_model::tri_3()");
			return n[i];
		}
	};

	std::vector <point_3> vp;
	std::vector <tri_3> vt;		//各个三角形旋转方向是朝外的，按照右手定则
	std::vector <int> vpli;		//为每个三角形分配一个整数
	std::vector <bool> vtri_ori;		//记录每个三角形的转向，其实可以通过vt计算得到,true:vt的外法线方向和平面vpli[i]一样

	stl_model () {};
	int init_stl (std::ifstream &infile);
	vector_3 get_normal (const int & n1, const int &n2, const int &n3) const;
	vector_3 get_normal (const tri_3 &tri) const {return get_normal (tri[0],tri[1],tri[2]);}
	vector_3 get_normal (const int &i) const {return get_normal(vt[i]);}
	plane_3 get_plane (const int &i) const {return plane_3 (vp[vt[i].n[0]],get_normal(i));}
	bbox_3 init_bbox (const int&i, const double &tol=0.000000001) const
	{
		using CGAL::to_double;
		FT minx (vp[vt[i].n[0]].x()), maxx (vp[vt[i].n[0]].x());
		FT miny (vp[vt[i].n[0]].y()), maxy (vp[vt[i].n[0]].y());
		FT minz (vp[vt[i].n[0]].z()), maxz (vp[vt[i].n[0]].z());
		for (int j(0);j!=3;++j)
		{
			if (vp[vt[i].n[j]].x() > maxx) maxx = vp[vt[i].n[j]].x();
			if (vp[vt[i].n[j]].x() < minx) minx = vp[vt[i].n[j]].x();
			if (vp[vt[i].n[j]].y() > maxy) maxy = vp[vt[i].n[j]].y();
			if (vp[vt[i].n[j]].y() < miny) miny = vp[vt[i].n[j]].y();
			if (vp[vt[i].n[j]].z() > maxz) maxz = vp[vt[i].n[j]].z();
			if (vp[vt[i].n[j]].z() < minz) minz = vp[vt[i].n[j]].z();
		}
		return bbox_3 (to_double(minx)-tol, to_double(miny) - tol, to_double(minz) -tol,
			to_double (maxx)+tol, to_double (maxy) + tol, to_double (maxz) + tol);
	}


	//检查是否有外露边界和非流形边界。并且检查各个三角面的旋转方向是否一致。未来还可以加上修复各个面朝向的部分，
	int check () const;	
	void clear () {vp.clear(); vt.clear(); vpli.clear ();}

	//计算该结构与给定平面的交集。所生成的loop记录在vl中，并且loop的旋转方向是view against pl的正向，用单面切割，并且只计算在模型内部的点
	//返回值是0：成功，否则弹出异常。如果没有交集，则vl和vp都是空的
	int plane_intersect (const plane_3 & pl, std::vector <point_2> &vp, std::vector<index_l> &vl) const;
	//vl中各个loop的方向为从pl的正向看下去的方向
	int plane_intersect (const plane_3 & pl, std::vector <point_3> &vp, std::vector <index_l>&vl) const;

	int plane_intersect (const plane_3 & pl, csection &sec)
	{
		return plane_intersect(pl,sec.vp,sec.vl);
	}
};



class frac_base
{
protected:
	frac_base ()
		: frac_id(-1),
		  plane_id(-1),
		  cohesion(0),
		  f_angle(0),
		  aperture(0){}

public:
	int		plane_id;			//裂隙所在平面的ID
	int		frac_id;			//裂隙的ID,正数
	FT		cohesion;
	FT		f_angle;
	FT		aperture;
};

class poly_frac: public frac_base
{
public:
	enum Type {Pos_side, Neg_side, Both_side};
	typedef std::vector <point_3> f_loop;
	std::vector <f_loop> vl;		//表示多边形裂隙的边界loop，vl[0]是外表面边界，其余是内表面边界，从plane_id这个平面看下去，外表面是逆时针，其余表面边界是顺时针
	Type type;						//记录这个裂隙的哪一侧参与切割，默认是Both_side

	void proj_plane (const plane_3& pl)		//将这个裂隙上的点都投影到pl上
	{
		for (int i(0);i!=(int)vl.size ();++i)
			for (int j(0);j!=(int)vl[i].size ();++j)
				vl[i][j] = pl.projection (vl[i][j]);
	}

	int to_2d (const plane_3&pl, std::vector<std::vector<point_2>> &vvp) const;
	int to_2d (const plane_3&pl, polygon_with_holes_2 &pwh) const;
	void revise_loops ()
	{
		for (int i(0);i!=(int)vl.size ();++i)
			std::reverse (vl[i].begin (),vl[i].end ());
	}
	bbox_3 init_bbox (double tol = 0.0000001) const
	{
		if (vl.empty ())
			throw std::logic_error ("init_bbox, poly is empty");
		using CGAL::to_double;
		FT minx (vl[0][0].x()), maxx (vl[0][0].x());
		FT miny (vl[0][0].y()), maxy (vl[0][0].y());
		FT minz (vl[0][0].z()), maxz (vl[0][0].z());
		for (int i(0);i!=vl.size ();++i)
			for (int j(0);j!=vl[i].size ();++j)
			{
				if (vl[i][j].x() > maxx) maxx = vl[i][j].x();
				if (vl[i][j].x() < minx) minx = vl[i][j].x();
				if (vl[i][j].y() > maxy) maxy = vl[i][j].y();
				if (vl[i][j].y() < miny) miny = vl[i][j].y();
				if (vl[i][j].z() > maxz) maxz = vl[i][j].z();
				if (vl[i][j].z() < minz) minz = vl[i][j].z();
			}
		return bbox_3 (to_double(minx)-tol, to_double(miny) - tol, to_double(minz) -tol,
			to_double (maxx)+tol, to_double (maxy) + tol, to_double (maxz) + tol);
	}
	void clear ()
	{ vl.clear ();}
};

class disc_frac: public frac_base
{
public:
	point_3 center;
	FT r;
	vector_3 normal;
};


template <typename K>
inline int add_set (std::vector<K> &p, const K &po)
{
	int i(0);
	for (i;i!=p.size ();++i)
		if (p[i] == po)
			break;
	if (i==p.size ())
		p.push_back (po);
	return i;
}

inline int add_set_pl (std::vector<plane_3> &vp, const plane_3 &p)
{
	int i(0);
	for (i;i!=vp.size ();++i)
		if (CGAL::parallel(vp[i],p) && vp[i].has_on(p.point()))
			break;

	if (i==vp.size ())
		vp.push_back (p);
	return i;
}

//在ve中添加新的线段，要求是如果e与ve中的其他线段相交， 则将这些相交的线段拆开。
int add_set_seg_3 (std::vector<segment_3> &ve, const segment_3 &e,std::vector <bbox_3> &vebox);

template <typename K>
inline void remove_set (std::vector<K> &vi, const std::vector <K> &tar)
{
	for (int i(0);i!=tar.size ();++i)
		vi.resize (std::remove (vi.begin (),vi.end (),tar[i]) - vi.begin ());
}

template <typename K>
inline void remove_set (std::vector<K> &vi, const K &tar)
{
	vi.resize (std::remove (vi.begin (),vi.end (),tar) - vi.begin ());
}

//将一个倾向倾角表示的方向计算出来，当然这种计算不是严格的，因为中间用到了sin和cos函数
inline vector_3 dip2normal (const FT& dip_dir_, const FT& dip_,const FT& dx_)
{
	double dip_dir = CGAL::to_double(dip_dir_*M_PI/180.0);
	double dip = CGAL::to_double(dip_*M_PI/180.0);
	double dx = CGAL::to_double(dx_*M_PI/180.0);
	return vector_3 (std::sin(dip)*std::cos(dip_dir-dx),
					-std::sin(dip)*std::sin(dip_dir-dx),
					std::cos(dip));
}

class right_most_compare_2 
{
	const std::vector<point_2> &vp;
	const edge &e;
public:
	right_most_compare_2 (const std::vector<point_2> & vp_, const edge&e_): vp(vp_), e(e_) {}
	bool operator() (const edge&e1, const edge&e2) const		//e1在e2右侧时返回true,即e1比e2对应的右角小
	{
		vector_2 v (vp[e.n[1]],vp[e.n[0]]);
		vector_2 v1 (vp[e1.n[0]],vp[e1.n[1]]);
		vector_2 v2 (vp[e2.n[0]],vp[e2.n[1]]);

		FT cs_prod1 (v.x()*v1.y()-v1.x()*v.y());
		FT cs_prod2 (v.x()*v2.y()-v2.x()*v.y());

		FT dot_prod1(v*v1);
		FT dot_prod2(v*v2);

		//    3    |    2
		//         |
		//--------------------
		//    4    |    1
		//         |   
		//         v
		int quad1, quad2;	//象限分布按照上面的顺序来,如果将v的方向视为y轴逆向
		if (cs_prod1>=0 && dot_prod1>0)
			quad1 = 1;
		else if (cs_prod1 >0 && dot_prod1<=0)
			quad1 = 2;
		else if (cs_prod1<=0 && dot_prod1<0)
			quad1 = 3;
		else
			quad1 = 4;

		if (cs_prod2>=0 && dot_prod2>0)
			quad2 = 1;
		else if (cs_prod2 >0 && dot_prod2<=0)
			quad2 = 2;
		else if (cs_prod2<=0 && dot_prod2<0)
			quad2 = 3;
		else
			quad2 = 4;

		if (quad2 > quad1)
			return true;		//e1在e2右侧
		
		if (quad1 > quad2)
			return false;		//e1在e2左侧

		FT scs_prod1 (CGAL::square(dot_prod1)*v2.squared_length ());
		FT scs_prod2 (CGAL::square(dot_prod2)*v1.squared_length ());

		if (quad1 == 1 || quad1 == 3)
			return (scs_prod1 > scs_prod2);		//返回true表明e1在e2右侧
		else
			return (scs_prod1 < scs_prod2);		//返回true表明e1在e2右侧
	}
};

class right_most_compare_3
{
	const std::vector<point_3> &vp;
	const plane_3 &pl;
	const edge &e;
public:
	right_most_compare_3 (const plane_3 &pl_, const std::vector<point_3> &vp_, const edge&e_)
		: pl(pl_), vp(vp_), e(e_) {}
	bool operator () (const edge &e1, const edge &e2) const		//从pl的正向往下看，e1在e2右侧时返回true,即e1比e2对应的右角小
	{

		vector_3 v (vp[e.n[1]], vp[e.n[0]]);
		vector_3 v1 (vp[e1.n[0]],vp[e1.n[1]]);
		vector_3 v2 (vp[e2.n[0]],vp[e2.n[1]]);

		FT cs_prod1 (CGAL::cross_product(v,v1) * pl.orthogonal_vector());
		FT cs_prod2 (CGAL::cross_product(v,v2) * pl.orthogonal_vector());

		FT dot_prod1 (v*v1);
		FT dot_prod2 (v*v2);

		int quad1,quad2;

		if (cs_prod1>=0 && dot_prod1>0)
			quad1 = 1;
		else if (cs_prod1 > 0 && dot_prod1<=0)
			quad1 = 2;
		else if (cs_prod1 <=0 && dot_prod1<0)
			quad1 = 3;
		else
			quad1 = 4;

		if (cs_prod2>=0 && dot_prod2>0)
			quad2 = 1;
		else if (cs_prod2 >0 && dot_prod2<=0)
			quad2 = 2;
		else if (cs_prod2<=0 && dot_prod2<0)
			quad2 = 3;
		else
			quad2 = 4;

		if (quad2 > quad1)
			return true;		//e1在e2右侧
		
		if (quad1 > quad2)
			return false;		//e1在e2左侧

		FT scs_prod1 (CGAL::square(dot_prod1)*v2.squared_length ());
		FT scs_prod2 (CGAL::square(dot_prod2)*v1.squared_length ());

		if (quad1 == 1 || quad1 == 3)
			return (scs_prod1 > scs_prod2);		//返回true表明e1在e2右侧
		else
			return (scs_prod1 < scs_prod2);		//返回true表明e1在e2右侧
	}
};

class right_most_compare_3_vector
{
	const vector_3 &n;
	const vector_3 &org;
public:
	right_most_compare_3_vector (const vector_3 &n_, const vector_3 &org_) : n(n_) , org(org_){}
	bool operator () (const vector_3 &v1, const vector_3 &v2) const		//逆着n看下去，v1在v2右侧时返回true，即v1比v2对应的右角小
	{
		vector_3 v (-org);
		FT cs_prod1 (CGAL::cross_product(v,v1) * n);
		FT cs_prod2 (CGAL::cross_product(v,v2) * n);

		FT dot_prod1 (v*v1);
		FT dot_prod2 (v*v2);

		int quad1,quad2;

		if (cs_prod1>=0 && dot_prod1>0)
			quad1 = 1;
		else if (cs_prod1 > 0 && dot_prod1<=0)
			quad1 = 2;
		else if (cs_prod1 <=0 && dot_prod1<0)
			quad1 = 3;
		else
			quad1 = 4;

		if (cs_prod2>=0 && dot_prod2>0)
			quad2 = 1;
		else if (cs_prod2 >0 && dot_prod2<=0)
			quad2 = 2;
		else if (cs_prod2<=0 && dot_prod2<0)
			quad2 = 3;
		else
			quad2 = 4;

		if (quad2 > quad1)
			return true;		//e1在e2右侧
		
		if (quad1 > quad2)
			return false;		//e1在e2左侧

		FT scs_prod1 (CGAL::square(dot_prod1)*v2.squared_length ());
		FT scs_prod2 (CGAL::square(dot_prod2)*v1.squared_length ());

		if (quad1 == 1 || quad1 == 3)
			return (scs_prod1 > scs_prod2);		//返回true表明e1在e2右侧
		else
			return (scs_prod1 < scs_prod2);		//返回true表明e1在e2右侧
	}
};

class edge_compare
{
public:
	bool operator () (const edge&e1, const edge&e2)
	{
		if (e1.n[0] < e2.n[0])
			return true;
		if (e1.n[0] == e2.n[0])
			if (e1.n[1] < e2.n[1])
				return true;
		return false;
	}
};

class point_compare_3 
{
public: 
	 bool operator () (const point_3&p1, const point_3&p2)
	 {
		 if (CGAL::compare_xyz(p1,p2) == CGAL::Comparison_result::SMALLER)
			 return true;
		 else
			 return false;
	 }
};
//计算np中索引所指向的点组成的法线方向。假设pn中所有点都在pl上（exactly）
vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np, const plane_3 &pl);
//同上，计算np中索引所指向的点组成的法线方向。假设pn中所有点都（exactly）在同一个平面上
vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np);
//同上，计算vp中所有点的组成边的法线方向。假设vp中所有点都（exactly）在同一个平面上
vector_3 oriented_normal (const std::vector<point_3> &vp);

//检查tar是否在向量beg，end的左侧，返回1表示在左侧，返回0表示在向量上，返回2表示在右侧
inline int is_left_2(const point_2 &beg, const point_2 &end, const point_2 &tar)
{
	FT temp = (end.x()-beg.x()) * (tar.y() - beg.y()) - (tar.x()-beg.x())*(end.y()-beg.y());
	if (temp>0)
		return 1;
	if (temp<0)
		return 2;
	else
		return 0;
}

//同上，假设beg和end都在pl上，结果是viewed against pl的正向
inline int is_left_3(const plane_3 &pl, const point_3 &beg, const point_3 &end, const point_3 &tar)
{
	FT temp = CGAL::cross_product (vector_3(tar,beg),vector_3(tar,end)) * pl.orthogonal_vector();
	if (temp > 0)
		return 1;
	if (temp < 0)
		return 2;
	else 
		return 0;
}

//返回0表示在多边形上，返回1表示在多边形内，返回2表示在多边形外
int point_in_polygon_2 (const std::vector<point_2> &vp, const std::vector<index_l> & vl, const point_2& tar);

////前提是ob和tar不在边界处相交，只在顶点处相交(不检查)。
////1表示tar的顶点完全包含在ob内，返回2表示tar有一个顶点包含在ob外
//int polygon_in_polygon_2 (const std::vector <point_2> &vp, const index_l &ob, const index_l &tar);
//同上，假设vl中的点都在pl上,并且从pl正向往下看逆时针的loop是内部，顺时针的loop是外部
//要求有且仅有一个外表面loop
//返回0表示在多边形上,返回1表示在多边形内，返回2表示在多边形外
int point_in_polygon_3 (const std::vector<point_3> &vp, const std::vector<index_l> & vl,const plane_3 &pl, const point_3 &tar);
int point_in_polygon_3 (const std::vector<point_3> &vp, const std::vector<loop> & vl,const plane_3 &pl, const point_3 &tar);
int point_in_polygon_3 (const std::vector<std::vector<point_3>> &vl, const plane_3&pl, const point_3 &tar);
int point_in_polygon_3 (const std::vector<point_3> &vp, const index_l & l,const plane_3 &pl, const point_3 &tar);

//判断一个点和一个多面体的关系，假设b中的各个面围成一个封闭空间（体积大于0），并且忽略各个面的方向性
//0：这个点在多面体上
//1：这个点在多面性内
//2：这个点在多边形外
int point_in_polyhedron_3 (const std::vector <point_3> &vp, const std::vector <face> &vfa, const std::vector<plane_3> &vpl, const block &b, const std::vector<bbox_3> &vbbox, const point_3 &tar);

//判断一个线段和一个多边形的位置关系
//-1:线段和多边形不相交
//1：线段和多边形共面(不管相交不相交)
//2：线段的某个端点在多边形的边界上
//3：线段的某个端点在多边形内
//4：线段和多边形的边界相交
//5：线段和多边形内相交
int seg_intersect_polygon_3 (const std::vector<point_3> &vp, const std::vector<loop> &vl, const plane_3 &pl, const point_3 &beg, const point_3 &end);

//-1:不相交
//1:可能相交
inline int seg_intersect_bbox_3 (const bbox_3 &bbox, const point_3 &beg, const point_3 &end)
{
	for (int i(0);i!=3;++i)
	{
		if (beg[i]<bbox.min(i) && end[i]<bbox.min(i))
			return -1;
		if (beg[i]>bbox.max(i) && end[i]>bbox.max(i))
			return -1;
	}
	return 1;
}

//前提是ob和tar不在边界处相交，只在顶点处相交(不检查)。
//要求ob是逆时针旋转
//1表示tar的顶点完全包含在ob内，返回2表示tar有一个顶点包含在ob外
int polygon_in_polygon_3 (const std::vector <point_3> &vp, const index_l & ob, const index_l &tar, const plane_3 &pl);

//前提是ob和tar不在边界处相交，只在顶点处相交,并且ob和tar除了边界处没有重合的地方
//要求ob和tar都围成一个封闭的空间。体积不为0，判断的时候忽略其方向性。并且需要给出已经计算出的体积
//返回1表示tar的顶点完全包含在ob内，返回2表示tar至少有一个顶点包含在ob外
int polyhedron_in_polyhedron_3 (const std::vector<point_3> &vp, 
								const std::vector<face> &vfa,
								const std::vector<plane_3> &vpl,
								const block&ob, const block& tar, 
								const FT &v_ob, const FT &v_tar, 
								const std::vector<bbox_3>& vbox_ob, const std::vector<bbox_3>& vbox_tar);

//假设vl上的点都在同一个平面pl内
point_3 get_point_in_face_3 (const std::vector<point_3> &vp, const plane_3 &pl, const face &fa);

//用户需要保证blk的体积不为零，并且blk围成一个封闭的空间，如果blk的体积小于零，则调用时需将reverse设置为true
point_3 get_point_in_polyhedron_3 (const std::vector<point_3> &vp, 
								   const std::vector <face> &vfa,
								   const std::vector <plane_3> &vpl, 
								   const block&blk, const std::vector<bbox_3> &vbbox,const bool& reverse=false);

//分裂一个loop，同时分裂这个loop中所包含的边界信息
void split_loop (const loop& l, std::vector<loop> &vl);

//用po分裂vseg中的线段，并且更新vbox_seg
int split_seg (std::vector<segment_3> &vseg, std::vector<bbox_3> &vbox_seg, const point_3 & po);

//为blk的每一个边界初始化一个边界盒
void get_vbox_of_blk (const std::vector<point_3> &vp, const std::vector<face> &vfa, const block &blk, std::vector<bbox_3> &vbox);

//inline void out_line_3_rhino (std::ostream &os, const plane_3& pl, const std::vector<point_2> &vp, const std::vector<loop> &vl)
//{
//	for (int i(0);i!=vl.size ();++i)
//	{
//		os<<"polyline ";
//		for (int j(0);j!=vl[i].vpi.size ();++j)
//		{
//			point_3 temp (pl.to_3d(vp[vl[i].vpi[j]]));
//			os<<temp.x()<<','<<temp.y()<<','<<temp.z()<<std::endl;
//		}
//		os<<"c"<<std::endl;
//	}
//}

inline void out_line_3_rhino (std::ostream &os, const std::vector<point_3> &vpo, const index_l&l)
{
	os<<"polyline ";
	for (int j(0);j!=l.size ();++j)
	{ 
		const point_3 &temp (vpo[l[j]]);
		os<<temp.x()<<','<<temp.y()<<','<<temp.z()<<std::endl;
	}
	os<<"c"<<std::endl;
}

inline void out_line_3_rhino (std::ostream &os, const segment_3 &s)
{
	os<<"line "
		<<s.source ().x()<<","<<s.source ().y()<<","<<s.source ().z()<<" "
		<<s.target ().x()<<","<<s.target ().y()<<","<<s.target ().z()<<" "<<std::endl ;
}

//输出一个平面的stl，flag=true，将该平面输出成一个单独的stl文件
void out_face_STL (std::ostream &outfile, const std::vector <point_3> &vp, const face &fa, const plane_3 &pl, bool flag ,const std::string &s = "");

//计算两个pf的交点，rese记录所有的线段，resp记录所有的点，返回值是（rese增加的值，resp增加的值）
//计算出来的交点不会是rese中线段的端点或在这些线段上
//添加rese时用的是add_set_seg_3,添加resp时是直接添加
std::pair<int,int> intersection (const poly_frac &pf1, const plane_3 &pl1,
								 const poly_frac &pf2, const plane_3 &pl2,
								 std::vector<segment_3> & rese, std::vector<point_3> &resp,
								 std::vector <bbox_3> &vebox);


//计算loop和直线(p, v)的交点，返回值保存在res中，表示计算出来的交点，即p+res[i]*v
//precondition: loop和给定直线共面
//如果loop的某个边界线段和直线相交，则返回两个边界线段的顶点
void intersection_points (const std::vector<point_3> & loop, const point_3&p, const vector_3& v, std::vector<FT> &res);
void intersection_points (const std::vector<std::vector<point_3>> &vl, const point_3&p, const vector_3& v, std::vector<FT> &res);
inline void intersection_points (const poly_frac &pf, const point_3&p, const vector_3&v, std::vector <FT>&res)
{ intersection_points(pf.vl,p,v,res);}

//判断segs给定的线段是否包含在pf中或在其边界上，res[i]中值表示线段segs[i],segs[i+1]中的线段是否在pf上.
//pf中各个loop的旋转方向是从pl的正向看下去的方向
//要求segs中的元素不重复
void mark_segments (const poly_frac &pf, const plane_3 &pl, const point_3&p, const vector_3&v, const std::vector <FT> &segs, std::vector <bool> &res);

//将返回值记为x，则有tar = p+x*v
inline FT get_para (const point_3&p, const vector_3&v, const point_3& tar)
{
	using CGAL::abs;
	if (abs(v.x()) >= abs(v.y()) && abs(v.x()) >= abs (v.z()))
		return (tar.x()-p.x())/v.x();
	else if (abs(v.y()) >= abs(v.x()) && abs(v.y()) >= abs(v.z()))
		return (tar.y()-p.y())/v.y();
	else
		return (tar.z()-p.z())/v.z();
}

inline bbox_3 seg2bbox (const segment_3&seg, const double &tol = 0.00000001)
{
	return bbox_3 (CGAL::to_double (CGAL::min(seg.source ().x(), seg.target ().x())) - tol,
		CGAL::to_double (CGAL::min(seg.source ().y(), seg.target ().y())) - tol,
		CGAL::to_double (CGAL::min(seg.source ().z(), seg.target ().z())) - tol,
		CGAL::to_double (CGAL::max(seg.source ().x(), seg.target ().x())) + tol,
		CGAL::to_double (CGAL::max(seg.source ().y(), seg.target ().y())) + tol,
		CGAL::to_double (CGAL::max(seg.source ().z(), seg.target ().z())) + tol);
}

inline bbox_3 point2bbox (const point_3&p, const double &tol=0.00000001)
{
	using CGAL::to_double;
	return bbox_3 (to_double (p.x()) - tol, to_double(p.y())-tol, to_double (p.z()) - tol,
		to_double (p.x()) + tol, to_double (p.y()) +tol, to_double (p.z()) + tol);
}

//将一个圆盘转换成一个n边形，当然，转换不可能是exact的
void disc2pf (const point_3& center, const FT &r, const plane_3 &pl, poly_frac & res, const int& n = 10);

//用来标记三角形面，在输出outer_boundary和inner_boundary的STL时有用到。只能用来标记嵌套过一次的面
template <typename CDT>
inline void mark_domains (CDT &cdt, typename CDT::Face_handle start, int index, std::list<typename CDT::Edge>& border)
{
	if (start->info().nesting_level !=-1)
		return ;
	std::list<CDT::Face_handle> queue;
	queue.push_back (start);
	while (!queue.empty())
	{
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1)
		{
			fh->info().nesting_level = index;
			for (int i(0);i<3;++i)
			{
				CDT::Edge e(fh,i);
				CDT::Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1)
					if (cdt.is_constrained(e)) border.push_back (e);
					else queue.push_back(n);
			}
		}
	}
}


//和上面那个函数结合使用，（其实是上面那个函数和这个函数结合使用）
template <typename CDT>
inline void mark_domains (CDT &cdt)
{
	for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it!=cdt.all_faces_end();++it)
		it->info().nesting_level = -1;
	std::list<CDT::Edge> border;
	mark_domains(cdt,cdt.infinite_face(),0,border);
	while (!border.empty())
	{
		CDT::Edge e = border.front();
		border.pop_front();
		CDT::Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1)
			mark_domains(cdt,n,e.first->info().nesting_level+1,border);
	}
}

//这个函数用来给一个Constrained Delaunay Triangulation (CDT)中添加多边形限制的
template <typename CDT>
inline void insert_polygon(CDT& cdt, const std::vector<point_3> &vp, const plane_3& pl, const loop& l)
{
	if (l.vpi.empty())
		return;
	CDT::Vertex_handle v_prev = cdt.insert (pl.to_2d (vp[l.vpi.back()]));
	for (int i(0);i!=l.vpi.size ();++i)
	{
		CDT::Vertex_handle vh = cdt.insert (pl.to_2d (vp[l.vpi[i]]));
		cdt.insert_constraint (vh,v_prev);
		v_prev = vh;
	}
}

}