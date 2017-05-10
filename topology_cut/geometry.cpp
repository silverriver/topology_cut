#include "geometry.h"
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <CGAL\Constrained_Delaunay_triangulation_2.h>
#include <CGAL\Triangulation_face_base_with_info_2.h>
using namespace std;

namespace TC
{

bool operator == (const index_l &l1, const index_l &l2)
{

	if (l1.size () != l2.size ())
		return false;
	for (int i(0);i!=l1.size ();++i)
	{
		bool flag (false);
		for (int j(0);j!=l2.size ();++j)
		{
			if (l1[i] == l2[j])			//如果找到两个相同的元素
			{
				flag = true;		//vpi[i]等于l.vpi中的某个元素
				int count(0);
				for (; count!=l1.size ();++count)
				{
					int ic (i+count), jc (j+count);
					if (ic>=l1.size ()) ic-=(int)l1.size ();
					if (jc>=l1.size ()) jc-=(int)l1.size ();
					if (l1[ic] != l2[jc])
						break;
				}
				if (count == l1.size ())
					return true;		//成功结束
			}
		}
		if (flag)
			return false;		//如果vpi[i]不等于l.vpi中的任何一个元素
	}
	return false;
}

int add_set_seg_3 (std::vector<segment_3> &ve, const segment_3 &e,std::vector <bbox_3> &vebox)
{
	if (e.is_degenerate ())
		return 0;

	if (ve.empty ())
	{
		ve.push_back (e);
		vebox.push_back (seg2bbox(e));
		return 1;
	}

	vector<FT> vep;
	vector<pair<FT,FT>> vpair;
	vep.reserve (10);
	vpair.reserve (10);
	vep.push_back (0);
	vep.push_back (1);
	vector_3 tempv (e.source (), e.target());
	FT occupied_length(0);		//记录已e上已经有多少已经被添加的部分
	size_t ve_size = ve.size ();
	bbox_3 box (seg2bbox(e));
	for (int i(0);i!=ve_size;++i)		//新添加的边不用判断
	{
		if (!CGAL::do_overlap (vebox[i],box))
			continue;

		CGAL::cpp11::result_of<K::Intersect_3(segment_3,segment_3)>::type
			res = CGAL::intersection (ve[i],e);
		if (res)
		{
			if (const point_3* p = boost::get<point_3> (&*res))		//相交的是一个点
			{
				vep.push_back (get_para(e.source(),tempv,*p));

				if ( ((*p) != ve[i].source()) && ((*p) != ve[i].target ()))		//e与ve[i]的内点相交,将ve[i]分成两个
				{
					ve.push_back (segment_3 (*p,ve[i].target ()));
					ve[i] = segment_3 (ve[i].source(),*p);
					vebox.push_back (seg2bbox(ve.back ()));
					vebox[i] = seg2bbox (ve[i]);
				}
			}
			if (const segment_3* new_e = boost::get<segment_3> (&*res))		//相交的是一个线段
			{
				vpair.push_back (pair<FT,FT> (get_para(e.source (),tempv,new_e->source ()),
											get_para (e.source (), tempv, new_e->target ())));
				occupied_length = occupied_length + CGAL::abs (vpair.back ().first - vpair.back ().second );	//因为ve中的边都没有重合的
				vep.push_back (vpair.back ().first );
				vep.push_back (vpair.back ().second );
				//将ve[i]分裂
				vector_3 ve_tempv (ve[i].source (),ve[i].target ());
				FT p1_ (get_para(ve[i].source(),ve_tempv,new_e->source ()));
				FT p2_ (get_para(ve[i].source(),ve_tempv,new_e->target ()));
				if (CGAL::abs(p1_-p2_) == 1)
					continue;		//ve与e重合
				FT p1 (CGAL::min (p1_,p2_)), p2 (CGAL::max (p1_,p2_));		//计算出来的新边端点在ve[i]上对应的参数
				
				if (p1 == 0)		//最多新生成s1, s2 和*new_e这几个线段
				{
					if (p2 == 1)
						continue;
					else
					{
						ve[i] = *new_e;
						ve.push_back (segment_3 (ve[i].source () + p2*ve_tempv, ve[i].target ()));
						vebox[i] = seg2bbox(ve[i]);
						vebox.push_back (seg2bbox(ve.back ()));
					}
				}
				else if (p2 == 1)
				{
					ve[i] = *new_e;
					ve.push_back (segment_3 (ve[i].source (), ve[i].source () + p1*ve_tempv));
					vebox[i] = seg2bbox(ve[i]);
					vebox.push_back (seg2bbox(ve.back ()));
				}
				else
				{
					ve[i] = *new_e;
					vebox[i] = seg2bbox(ve[i]);
					ve.push_back (segment_3 (ve[i].source (), ve[i].source () + p1*ve_tempv));
					vebox.push_back (seg2bbox(ve.back ()));
					ve.push_back (segment_3 (ve[i].source () + p2*ve_tempv, ve[i].target ()));
					vebox.push_back (seg2bbox(ve.back ()));
				}

				//上面这部分代码可能更快一点
				//segment_3 s1 (ve[i].source (), ve[i].source () + p1*ve_tempv), s2 (ve[i].source () + p2*ve_tempv, ve[i].target ());		//最多新生成s1, s2 和*new_e这几个线段
				//ve[i] = *new_e;
				//vebox[i] = seg2bbox(ve[i]);
				//if (!s1.is_degenerate ())
				//{
				//	ve.push_back (s1);
				//	vebox.push_back (seg2bbox(ve.back ()));
				//}
				//if (!s2.is_degenerate ())
				//{
				//	ve.push_back (s2);
				//	vebox.push_back (seg2bbox(ve.back ()));
				//}
			}
		}
		if (occupied_length == 1)
			break;		//如果e上所有部分都已经被添加了，则直接跳出
	}

	if (occupied_length == 1)
		return (int) (ve.size ()-ve_size);		//不用添加e
	sort(vep.begin (), vep.end ());
	vep.resize (unique(vep.begin (), vep.end ()) - vep.begin ());
	if (vep.size () == 2)		//vep中只有两个端点，意思是ve中没有原始与e相交，直接添加
	{
		ve.push_back (e);
		vebox.push_back (seg2bbox(ve.back ()));
	}
	else
	{
		for (int i(0);i!=vep.size ();++i)
		{
			int i1(i+1);
			if (i1 >= vep.size ()) break;
			FT mid ((vep[i]+vep[i1])/2);
			int count(0);
			for (;count!=vpair.size ();++count)
				if (mid >= CGAL::min (vpair[count].first ,vpair[count].second ) && mid <= CGAL::max (vpair[count].first ,vpair[count].second ))
					break;
			if (count == vpair.size ())		//如果这一段还没有包含在ve中
			{
				ve.push_back (segment_3 (e.source () + vep[i]*tempv, e.source() + vep[i1]*tempv));
				vebox.push_back (seg2bbox(ve.back ()));
			}
		}
	}
	return (int) (ve.size ()-ve_size);
}

point_3 get_point_in_face_3 (const std::vector<point_3> &vp, const plane_3 &pl, const face &fa)
{
	if (fa.vl[0].vpi.size () <=2)
		throw logic_error ("get_point_in_face_3(), can not find interior points for degenerated points");

	for (int i(0);i!=fa.vl[0].vpi.size ();++i)
	{
		int i1 (i+1);
		if (i1>=fa.vl[0].vpi.size ()) i1 = 0;
		int i2 (i1+1);
		if (i2>=fa.vl[0].vpi.size ()) i2 = 0;
		point_3 mid (CGAL::midpoint(vp[fa.vl[0].vpi[i1]], vp[fa.vl[0].vpi[i2]]));
		if (point_in_polygon_3(vp,fa.vl,pl,CGAL::midpoint (vp[fa.vl[0].vpi[i]],mid)) == 1)
			return CGAL::midpoint (vp[fa.vl[0].vpi[i]],mid);

		line_3 temp_l (vp[fa.vl[0].vpi[i]],mid);
		vector<point_3> inter_point;		//记录temp_l和各个边界的交点
		inter_point.reserve (fa.vl.size ()*5);
		for (int j(0);j!=fa.vl.size ();++j)
		{
			for (int k(0);k!=fa.vl[j].vpi.size ();++k)
			{
				int k1(k+1);
				if (k1>=fa.vl[j].vpi.size ()) k1=0;
				const point_3 &p1 (vp[fa.vl[j].vpi[k]]), &p2 (vp[fa.vl[j].vpi[k1]]);
				CGAL::cpp11::result_of<K::Intersect_3(segment_3,line_3)>::type
					res = CGAL::intersection (segment_3 (p1,p2),temp_l);
				if (res)
				{
					if (const point_3*p = boost::get<point_3> (&*res))
						inter_point.push_back (*p);
				}
			}
		}
		point_compare_3 pcom;
		sort (inter_point.begin (), inter_point.end (), pcom);
		inter_point.resize (unique(inter_point.begin (),inter_point.end ()) - inter_point.begin ());

		if (inter_point.size () <=1)
			throw logic_error ("get_point_in_face_3, less than 2 intersection points");
		for (int j(0);j!=inter_point.size ();++j)
		{
			int j1(j+1);
			if (j1>=inter_point.size ()) break;
			if (point_in_polygon_3(vp,fa.vl,pl,CGAL::midpoint (inter_point[j],inter_point[j1])) == 1)
				return CGAL::midpoint (inter_point[j],inter_point[j1]);
		}
	}

	throw runtime_error ("get_point_in_face_3, bad luck~ no point generated");
	return point_3();
}

point_3 get_point_in_polyhedron_3 (const std::vector<point_3> &vp, 
								   const std::vector <face> &vfa, 
								   const std::vector <plane_3> &vpl, 
								   const block&blk, const vector<bbox_3> &vbbox,
								   const bool& reverse)
{
	double tol(1);
	for (int i(0);i!=vbbox.size ();++i)
		for (int j(0); j!=3;++j)
			if ((vbbox[i].max (j) - vbbox[i].min(j))<tol) tol=vbbox[i].max (j) - vbbox[i].min(j);

	FT r = pow (CGAL::to_double(tol), 1/3.0); 
	int count = 1;
	vector<point_3> vp_temp (blk.vfi.size ());
	vector<vector_3> vv_temp (blk.vfi.size ());
	while (count < 20)
	{
		for (int i(0);i!=blk.vfi.size ();++i)
		{
			if (count == 1)
			{
				vp_temp[i] = get_point_in_face_3(vp,vpl[vfa[blk.vfi[i].first].pli],vfa[blk.vfi[i].first]);
				vv_temp[i] = -vpl[vfa[blk.vfi[i].first].pli].orthogonal_vector();
				if (!(blk.vfi[i].second))
					vv_temp[i] = -vv_temp[i];
				if (!reverse)
					vv_temp[i] = -vv_temp[i];
				double r = sqrt (to_double(vv_temp[i].x()*vv_temp[i].x() + vv_temp[i].y()*vv_temp[i].y() + vv_temp[i].z()*vv_temp[i].z()));
				vv_temp[i] = (1/r)*vv_temp[i];
			}
			if (point_in_polyhedron_3(vp,vfa,vpl,blk,vbbox,vp_temp[i]+r*vv_temp[i]) == 1)
				return vp_temp[i]+r*vv_temp[i];
		}
		r = r/2.0;
	}
	throw logic_error ("get_point_in_polyhedron_3, bad luck!");
}

vector_3 stl_model::get_normal(const int & n1, const int &n2, const int &n3)const
{
	return CGAL::cross_product(vector_3(vp[n1],vp[n2]),vector_3 (vp[n1],vp[n3]));
}

int stl_model::init_stl (std::ifstream &infile)
{
	if (!infile)
		return 1;
	vp.clear(); vt.clear();
	string keyword1, keyword2;
	char temp('a');
	while (temp!='\n')
		infile.get (temp);

	while (infile)
	{
		infile>>keyword1>>keyword2;
		if (keyword1 == "endsolid")
			break;
		if (keyword1 != "facet")
			return 2;
		
		getline (infile, keyword1);		//忽略掉法向量
		infile>>keyword1>>keyword2;
		if (keyword1 !="outer" ||keyword2 !="loop")
			return 2;
		
		double x[3],y[3],z[3];
		for (int i(0);i!=3;++i)
		{
			infile>>keyword1;
			if (keyword1 !="vertex")
				return 2;
			infile>>x[i]>>y[i]>>z[i];
		}
		infile >>keyword1>>keyword2;
		if (keyword1 !="endloop" || keyword2!="endfacet")
			return 2;
		
		tri_3 tri;
		tri[0] = add_set (vp,point_3(x[0],y[0],z[0]));
		tri[1] = add_set (vp,point_3(x[1],y[1],z[1]));
		tri[2] = add_set (vp,point_3(x[2],y[2],z[2]));
		vt.push_back (tri);
	}
	return 0;
}

int stl_model::check () const
{
	vector<edge> ve;
	vector<int> ve_info;		//记录这个边缘是否是一个外露边缘
	vector<int> ve_info1;		//记录这个边缘是否是一个非流行边缘
	for (int i(0);i!=(int)vt.size ();++i)
	{
		if (vt[i][0] == vt[i][1] || vt[i][0] == vt[i][2] || vt[i][1] == vt[i][2])
			return 1;		//有重复顶点
		int count = (int)ve.size ();
		int n1 = add_set(ve,edge(vt[i][0],vt[i][1]));
		if (ve.size () > count)		//新添加了边
		{
			ve_info.push_back (0);
			ve_info1.push_back (0);
			count = (int)ve.size();
		}

		int n2 = add_set(ve,edge(vt[i][1],vt[i][2]));
		if (ve.size () > count)		//新添加了边
		{
			ve_info.push_back (0);
			ve_info1.push_back (0);
			count = (int)ve.size();
		}
		int n3 = add_set(ve,edge(vt[i][2],vt[i][0]));
		if (ve.size () > count)		//新添加了边
		{
			ve_info.push_back (0);
			ve_info1.push_back (0);
			count = (int)ve.size();
		}
		
		ve_info1[n1] ++;
		ve_info1[n2] ++;
		ve_info1[n3] ++;

		if (ve[n1].n[0] == vt[i][0])	//如果ve[n1]和edge(n0,n1)顺序一样
			ve_info[n1] ++;
		else
			ve_info[n1] --;

		if (ve[n2].n[0] == vt[i][1])	//如果ve[n2]和edge(n1,n2)顺序一样
			ve_info[n2] ++;
		else
			ve_info[n2] --;

		if (ve[n3].n[0] == vt[i][2])	//如果ve[n3]和edge(n2,n0)顺序一样
			ve_info[n3] ++;
		else
			ve_info[n3] --;
	}
	for (int i(0);i!=(int)ve.size ();++i)
	{
		if (ve_info[i] != 0 )
			return 2;
		if (ve_info1[i] != 2)
			return 3;
	}
	return 0;
}

vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np, const plane_3 &pl)
{
	if (np.size ()<3)
		throw logic_error("oriented_normal, point count less than 3");
	point_3 ori(pl.point());
	vector_3 res(0,0,0);
	for (int i(0);i!=np.size ();++i)
	{
		 int i1(i+1);
		 if (i1>=(int)np.size ())
			 i1 = 0;
		 res = res + CGAL::cross_product(vector_3(ori,vp[np[i]]), vector_3(ori,vp[np[i1]]));
	}
	return res;
}

vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np)
{
	if (np.size ()<3)
		throw logic_error("oriented_normal, point count less than 3");
	const point_3 &ori =vp[np[0]];
	vector_3 res(0,0,0);
	for (int i(1);i!=np.size ()-1;++i)
	{
		 int i1(i+1);
		 res = res + CGAL::cross_product(vector_3(ori,vp[np[i]]), vector_3(ori,vp[np[i1]]));
	}
	return res;
}

vector_3 oriented_normal (const std::vector<point_3> &vp)
{
	if(vp.size ()<3)
		throw logic_error ("oriented_normal, point count less than 3");
	const point_3 &ori = vp[0];
	vector_3 res(0,0,0);
	for (int i(1);i!=vp.size ()-1;++i)
	{
		 int i1(i+1);
		 res = res + CGAL::cross_product(vector_3(ori,vp[i]), vector_3(ori,vp[i1]));
	}
	return res;
}

int point_in_polygon_2 (const std::vector<point_2> &vp, const std::vector<index_l> & vl, const point_2& tar)
{
	//http://geomalgorithms.com/a03-_inclusion.html
	int wn(0);
	
	for (int i(0);i!=vl.size ();++i)
	{
		for (int j(0);j!=vl[i].size ();++j)
		{
			int j1 (j+1);
			if (j1 >= vl[i].size ())
				j1=0;
			int p0 = vl[i][j], p1 = vl[i][j1];
			if (segment_2 (vp[p0],vp[p1]).has_on (tar))		//在线段上
				return 0;
			if (vp[p0].y() <= tar.y())		
			{
				if (vp[p1].y() > tar.y())		//an upward crossing
					if (is_left_2 (vp[p0], vp[p1],tar) ==1)		//tar left of edge
						++wn;
			}
			else
			{
				if (vp[p1].y() <= tar.y())		//an downward crossing
					if (is_left_2 (vp[p0],vp[p1],tar) == 2)			//tar right of edeg
						--wn;
			}
		}
	}
	if (wn == 0)
		return 2;
	else
		return 1;
}

int point_in_polygon_3 (const std::vector<point_3> &vp, const std::vector<index_l> & vl,const plane_3 &pl, const point_3 &tar)
{
	int wn(0);
	plane_3 temp (tar,pl.base1 ());		//temp ⊥ pl
	for (int i(0);i!=vl.size ();++i)
	{
		for (int j(0);j!=vl[i].size ();++j)
		{
			int j1 (j+1);
			if (j1 >= vl[i].size ())
				j1 = 0;
			int p0 = vl[i][j], p1 = vl[i][j1];
			if (segment_3(vp[p0], vp[p1]).has_on(tar))
				return 0;
			CGAL::Oriented_side p0_side (temp.oriented_side(vp[p0])), p1_side (temp.oriented_side(vp[p1]));
			if (p0_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)
			{
				if (p1_side == CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an upward crossing
					if (is_left_3 (pl,vp[p0],vp[p1],tar) == 1)				//tar is left of edge
						++wn;
			}
			else
			{
				if (p1_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an downward crossing
					if (is_left_3 (pl,vp[p0],vp[p1],tar) == 2)
						--wn;
			}
		}
	}
	if (wn == 0)
		return 2;
	else
		return 1;

}

int point_in_polygon_3 (const std::vector<point_3> &vp, const std::vector<loop> & vl,const plane_3 &pl, const point_3 &tar)
{
	int wn(0);
	plane_3 temp (tar,pl.base1 ());		//temp ⊥ pl
	for (int i(0);i!=vl.size ();++i)
	{
		for (int j(0);j!=vl[i].vpi.size ();++j)
		{
			int j1 (j+1);
			if (j1 >= vl[i].vpi.size ())
				j1 = 0;
			int p0 = vl[i].vpi[j], p1 = vl[i].vpi[j1];
			if (segment_3(vp[p0], vp[p1]).has_on(tar))
				return 0;
			CGAL::Oriented_side p0_side (temp.oriented_side(vp[p0])), p1_side (temp.oriented_side(vp[p1]));
			if (p0_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)
			{
				if (p1_side == CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an upward crossing
					if (is_left_3 (pl,vp[p0],vp[p1],tar) == 1)				//tar is left of edge
						++wn;
			}
			else
			{
				if (p1_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an downward crossing
					if (is_left_3 (pl,vp[p0],vp[p1],tar) == 2)
						--wn;
			}
		}
	}
	if (wn == 0)
		return 2;
	else
		return 1;

}

int point_in_polygon_3 (const std::vector<std::vector<point_3>> &vl, const plane_3&pl, const point_3 &tar)
{
	int wn(0);
	plane_3 temp (tar, pl.base1 ());
	for (int i(0);i!=vl.size ();++i)
	{
		for (int j(0);j!=vl[i].size ();++j)
		{
			int j1 (j+1);
			if (j1 >= vl[i].size ())
				j1 = 0;
			const point_3 & p0 = vl[i][j], p1 = vl[i][j1];
			if (segment_3(p0, p1).has_on(tar))
				return 0;
			CGAL::Oriented_side p0_side (temp.oriented_side(p0)), p1_side (temp.oriented_side(p1));
			if (p0_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)
			{
				if (p1_side == CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an upward crossing
					if (is_left_3 (pl,p0,p1,tar) == 1)				//tar is left of edge
						++wn;
			}
			else
			{
				if (p1_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an downward crossing
					if (is_left_3 (pl,p0,p1,tar) == 2)
						--wn;
			}
		}
	}
	if (wn == 0)
		return 2;
	else
		return 1;
}

int point_in_polygon_3 (const std::vector<point_3> &vp, const index_l & l,const plane_3 &pl, const point_3 &tar)
{
	int wn(0);
	plane_3 temp (tar,pl.base1 ());		//temp ⊥ pl

	for (int j(0);j!=l.size ();++j)
	{
		int j1 (j+1);
		if (j1 >= l.size ())
			j1 = 0;
		int p0 = l[j], p1 = l[j1];
		if (segment_3(vp[p0], vp[p1]).has_on(tar))
			return 0;
		CGAL::Oriented_side p0_side (temp.oriented_side(vp[p0])), p1_side (temp.oriented_side(vp[p1]));
		if (p0_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)
		{
			if (p1_side == CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an upward crossing
				if (is_left_3 (pl,vp[p0],vp[p1],tar) == 1)				//tar is left of edge
					++wn;
		}
		else
		{
			if (p1_side != CGAL::Oriented_side::ON_POSITIVE_SIDE)		//an downward crossing
				if (is_left_3 (pl,vp[p0],vp[p1],tar) == 2)
					--wn;
		}
	}

	if (wn == 0)
		return 2;
	else
		return 1;
}

int point_in_polyhedron_3 (const std::vector <point_3> &vp, const std::vector <face> &vfa, const std::vector<plane_3> &vpl, const block &b, const std::vector<bbox_3> &vbbox, const point_3 &tar)
{
	using CGAL::to_double;
	double tol (0.00000001);
	double minx(vbbox[0].xmin()-tol),miny(vbbox[0].ymin()-tol),minz(vbbox[0].zmin()-tol);
	double maxx(vbbox[0].xmax()+tol),maxy(vbbox[0].ymax()+tol),maxz(vbbox[0].zmax()+tol), r;
	
	for (int i(0);i!=vbbox.size ();++i)
	{
		if (vbbox[i].xmax () > maxx) maxx = vbbox[i].xmax();
		if (vbbox[i].xmin () < minx) minx = vbbox[i].xmin();
		if (vbbox[i].ymax () > maxy) maxy = vbbox[i].ymax();
		if (vbbox[i].ymin () < miny) miny = vbbox[i].ymin();
		if (vbbox[i].zmax () > maxz) maxz = vbbox[i].zmax();
		if (vbbox[i].zmin () < minz) minz = vbbox[i].zmin();
	}

	if(tar.x() > maxx || tar.x() < minx
		||tar.y () > maxy || tar.y() < miny
		||tar.z () > maxz || tar.z() < minz)
		return 2;

	r = sqrt ((maxx-minx)*(maxx-minx) + (maxy-miny)*(maxy-miny) + (maxz-minz)*(maxz-minz)) *1.1;
	srand (time(0));
	while (1)
	{
		int crossing = 0;
		double alpha (-M_PI/2.0 + rand()/double (RAND_MAX)*M_PI);
		double theta (rand()/double (RAND_MAX)*2*M_PI);
		point_3 tar_end(tar.x()+r*cos(alpha)*cos(theta), tar.y()+r*cos(alpha)*sin(theta), tar.z()+r*sin(theta));
		int count (0);
		for (; count !=b.vfi.size(); ++count);
		{
			if (seg_intersect_bbox_3(vbbox[count],tar,tar_end) == -1)
				continue;
			int code = seg_intersect_polygon_3(vp,vfa[b.vfi[count].first].vl,vpl[vfa[b.vfi[count].first].pli],tar, tar_end);
			if (code == -1)
				continue;
			else if (code == 1 || code == 4)	//如果这个ray和点或者边相交，则作为退化情况
				break;
			else if (code == 2 || code == 3)
				return 0;
			else if (code == 5)
				crossing++;
			else
				throw logic_error ("illegal return value from seg_intersect_polygon_3");
		}
		if (count != b.vfi.size ())
			continue;
		if ( (crossing %2) == 1)
			return 1;
		else 
			return 2;
	}
}

int seg_intersect_polygon_3 (const std::vector<point_3> &vp, const std::vector<loop> &vl, const plane_3 &pl, const point_3 &beg, const point_3 &end)
{
	FT nom = -pl.d () - pl.orthogonal_vector()*vector_3(CGAL::ORIGIN,beg);
	FT denom = pl.orthogonal_vector() * vector_3(beg,end);
	FT t; //beg + t(end - beg);
	if (denom == 0)
		if (nom == 0)
			return 1;		//共面
		else
			return -1;		//不相交
	else
		t = nom/denom;

	if (t < 0 || t > 1)
		return -1;		//不相交

	point_3 p_temp (beg+ t*vector_3(beg,end));
	int code = point_in_polygon_3(vp,vl,pl,p_temp);
	switch (code)
	{
	case 0:
		if (t == 0 || t == 1)
			return 2;
		else
			return 4;
		break;
	case 1:
		if (t == 0 || t == 1)
			return 3;
		else
			return 5;
		break;
	default:
		return -1;
	}
}

int polygon_in_polygon_3 (const std::vector <point_3> &vp, const index_l & ob, const index_l &tar, const plane_3 &pl)
{
	//判断tar是否完全包含在ob中
	//计算tar - ob是否为空,如果为空，则表明tar完全包含在ob中
	for (int i(0);i!=tar.size ();++i)
	{
		int i1(i+1);
		if (i1>=tar.size ())
			i1 = 0;
		if (point_in_polygon_3 (vp,ob,pl,CGAL::midpoint(vp[tar[i]], vp[tar[i1]])) == 2)		//tar的这个边在ob外
			return 2;		//表明tar没有完全包含在ob中
	}
	for (int i(0);i!=ob.size ();++i)
	{
		int i1(i+1);
		if (i1>=ob.size ())
			i1 = 0;
		if (point_in_polygon_3(vp,tar,pl,CGAL::midpoint(vp[ob[i]], vp[ob[i1]]))== 1)	//ob这个边包含在tar内	
			return 2;		//表明tar没有完全包含在ob中
	}
	return 1;

	//TODO：或者可以改进一下：只需要各取ob和tar的一个内点，判断这两个内点哪个包含在另一个多边形内，并且比较一下面积即可
}

int polyhedron_in_polyhedron_3 (const std::vector<point_3> &vp, const std::vector<face> &vfa, const std::vector<plane_3> &vpl, const block&ob, const block& tar, const FT& v_ob, const FT& v_tar, const vector<bbox_3>& vbox_ob, const vector<bbox_3>& vbox_tar)
{
	if (CGAL::abs (v_ob) < CGAL::abs (v_tar))
		return 2;
	if (v_ob == 0 || v_tar == 0)
		throw logic_error ("polyhedron_in_polyhedron_3, polyhedron have no volume");
	point_3 tar_p;
	if (v_ob > 0)
		tar_p = get_point_in_polyhedron_3(vp,vfa,vpl,tar,vbox_tar);
	else
		tar_p = get_point_in_polyhedron_3(vp,vfa,vpl,tar,vbox_tar,true);
	switch (point_in_polyhedron_3(vp,vfa,vpl,ob,vbox_ob,tar_p))
	{
	case 1:
		return 1;
	case 2:
		return 2;
	default:
		throw logic_error ("polyhedron_in_polyhedron_3, point on boundary");
	} 
}

void get_vbox_of_blk (const std::vector<point_3> &vp, const std::vector<face> &vfa, const block &blk, std::vector<bbox_3> &vbox)
{
	using CGAL::to_double;
	vbox.resize (blk.vfi.size ());
	double tol = 0.00000001;
	for (int i(0);i!=blk.vfi.size ();++i)
	{
		const face& f = vfa[blk.vfi[i].first];
		double minx,miny,minz,maxx,maxy,maxz;
		minx = maxx = to_double (vp[f.vl[0].vpi[0]].x());
		miny = maxy = to_double (vp[f.vl[0].vpi[0]].y());
		minz = maxz = to_double (vp[f.vl[0].vpi[0]].z());
		for (int j(0);j!=f.vl.size ();++j)
			for (int k(0);k!=f.vl[j].vpi.size ();++k)
			{
				if (minx < vp[f.vl[j].vpi[k]].x()) minx = to_double(vp[f.vl[j].vpi[k]].x());
				if (maxx > vp[f.vl[j].vpi[k]].x()) maxx = to_double(vp[f.vl[j].vpi[k]].x());
				if (miny < vp[f.vl[j].vpi[k]].y()) miny = to_double(vp[f.vl[j].vpi[k]].y());
				if (maxy > vp[f.vl[j].vpi[k]].y()) maxy = to_double(vp[f.vl[j].vpi[k]].y());
				if (minz < vp[f.vl[j].vpi[k]].z()) minz = to_double(vp[f.vl[j].vpi[k]].z());
				if (maxz > vp[f.vl[j].vpi[k]].z()) maxz = to_double(vp[f.vl[j].vpi[k]].z());
			}
			vbox[i] = bbox_3(minx - tol, miny - tol, minz - tol, maxx + tol, maxy + tol, maxz +tol);
	}
}

//int polygon_in_polygon_2 (const std::vector <point_2> &vp, const index_l &ob, const index_l &tar)
//{
//	vector<loop> temp;
//	temp.push_back (ob);
//	for (int i(0);i!=tar.vpi.size ();++i)
//		if (point_in_polygon_2(vp,temp,vp[tar.vpi[i]]) == 2)
//			return 2;
//	return 1;
//}

int poly_frac::to_2d (const plane_3&pl, std::vector<std::vector<point_2>> &vvp) const
{
	vvp.resize (vl.size ());	

	for (int i(0);i!=vl.size ();++i)
	{
		vvp[i].resize (vl[i].size ());		//这个变换不改变loop的旋转方向
		for (int j(0);j!=vl[i].size ();++j)
			vvp[i][j] = pl.to_2d(vl[i][j]);
	}
	return 0;
}

int poly_frac::to_2d (const plane_3&pl, polygon_with_holes_2 &pwh) const
{
	vector<vector<point_2>> vvp;
	to_2d(pl,vvp);
	pwh = polygon_with_holes_2(polygon_2(vvp[0].begin (),vvp[0].end ()));
	for (int i(1);i<vvp.size ();++i)
	{
		polygon_2 temp_p(vvp[i].begin (),vvp[i].end ());
		temp_p.reverse_orientation();
		pwh.add_hole (temp_p);
	}
	return 0;
}

int stl_model::plane_intersect (const plane_3 & pl, std::vector <point_2> &vp_2, std::vector<index_l> &vl) const
{
	vp_2.clear();
	vl.clear();
	vector<edge> ve;			//只记录在pl 正向上的三角形与pl的交线段
	vector<int> vet_same;		//用来记录每个边出现了几次与ve[i]同向的情况
	vector<int> vet_revise;		//用来记录每个边出现了几次与ve[i]逆向的情况

	for (int i(0);i!=vt.size ();++i)
	{
		int n_on(0), n_pos(0), n_neg(0);
		switch (pl.oriented_side(vp[vt[i].n[0]]))
		{
		case CGAL::ON_ORIENTED_BOUNDARY: n_on++; break;
		case CGAL::ON_POSITIVE_SIDE: n_pos++; break;
		case CGAL::ON_NEGATIVE_SIDE: n_neg++; break;
		}
		switch (pl.oriented_side(vp[vt[i].n[1]]))
		{
		case CGAL::ON_ORIENTED_BOUNDARY: n_on++; break;
		case CGAL::ON_POSITIVE_SIDE: n_pos++; break;
		case CGAL::ON_NEGATIVE_SIDE: n_neg++; break;
		}
		switch (pl.oriented_side(vp[vt[i].n[2]]))
		{
		case CGAL::ON_ORIENTED_BOUNDARY: n_on++; break;
		case CGAL::ON_POSITIVE_SIDE: n_pos++; break;
		case CGAL::ON_NEGATIVE_SIDE: n_neg++; break;
		}

		if (n_pos == 0)
			continue;		//如果三角形没有在pl正侧的点，则不参与计算
		if (n_on==3 || n_pos==3)
			continue;		//如果整个三角形都在pl正侧，或在pl中，则也不参与计算

		triangle_3 tri(vp[vt[i].n[0]], vp[vt[i].n[1]], vp[vt[i].n[2]]);
		CGAL::cpp11::result_of <K::Intersect_3(K::Plane_3,K::Triangle_3)>::type
			result = CGAL::intersection(pl,tri);
		if (!result)
			continue;
		if (const K::Segment_3 *s = boost::get<K::Segment_3> (&*result))	//交点是一个线段,交点是一个点的时候直接忽略
		{
			int beg = add_set(vp_2,pl.to_2d(s->source()));
			int end = add_set(vp_2,pl.to_2d(s->end ()));
			edge e;
			vector_3 v_temp (CGAL::cross_product(vector_3(*s),pl.orthogonal_vector()));
			if (v_temp * get_normal(i) <0)		//保证e的左侧是domain的内部
			{
				e.n[0] = end; e.n[1] = beg;
			}
			else
			{
				e.n[0] = beg; e.n[1] = end;
			}
			
			int ne = add_set (ve,e);
			if (vet_same.size () < ve.size ())		//如果这个edge是第一次被添加
			{
				vet_same.push_back (0);
				vet_revise.push_back (0);
			}

			if (ve[ne].n[0] == e.n[0])
				vet_same[ne]++;		//又出现了一次和ve[ne]同向的边
			else
				vet_revise[ne]++;	//又出现了一次和ve[ne]逆向的边
		}
	}
	for (int i(0);i!=ve.size ();++i)
	{
		if (vet_same[i] == 1 && vet_revise[i] == 0)		//normal situation
			continue;
		if (vet_same[i] != 1)					//因为假设stl没有非流行边界，所以每个边顶多出现两次，正反最多一次
			throw logic_error ("plane_intersect(), wrong stl format. may be nonmanifold edges");			
	}
	
	//如果某个边出现了正反各出现一次，那么其不参与生成loop
	int last_p(-1);		//记录目前loop最后一个edge，其实可以用一个bool类型替代。
	int count(0);
	for (int i(0);i!=ve.size ();++i)
		if (vet_revise[i] == 1)
			count++;
	while(1)
	{
		if (last_p == -1)
		{
			if (count == ve.size ())
				break;		//每个边都被添加了，皆大欢喜
			else
			{
				for (int i(0);i!=ve.size ();++i)		//寻找一个新的开始
					if (vet_revise[i] != 1)
					{
						vl.push_back(index_l());
						vl.back ().push_back (ve[i].n[0]);
						vl.back ().push_back (ve[i].n[1]);
						last_p = i;
						vet_revise[i] = 1;
						count++;
						break;
					}
			}
		}
		if (count == ve.size ())
			throw logic_error ("plane_intersection, loop is not finished");		//loop还没生成完成，边就没了

		vector<edge> proten_e;
		vector<int> e_index;
		for (int i(0);i!=ve.size ();++i)
		{
			if (vet_revise[i] == 1)
				continue;
			if (ve[i].n[0] == ve[last_p].n[1])
			{
				proten_e.push_back (ve[i]);
				e_index.push_back (i);
			}
		}
		if (proten_e.empty ())
			throw logic_error ("plane_intersection, can't find next edge");

		right_most_compare_2 temp_comp (vp_2,ve[last_p]);
		int right_most_n (0);
		for (int i(0);i!=proten_e.size ();++i)
			if (temp_comp (proten_e[right_most_n],proten_e[i]))
				right_most_n = i;
		
		last_p = e_index[right_most_n];
		vet_revise[last_p] = 1;
		count++;

		if (ve[last_p].n[1] == vl.back ()[0])
			last_p = -1;
		else
			vl.back ().push_back (ve[last_p].n[1]);
	}

	//vector<index_l> vl_temp;
	//vl_temp.reserve (vl.size ()+2);
	//for (int i(0);i!=vl.size ();++i)
	//	split_index_l(vl[i],vl_temp);
	//vl.swap (vl_temp);		//将每个loop都变成简单的loop
	return 0;
}

int stl_model::plane_intersect (const plane_3 & pl, std::vector <point_3> &vp_3, std::vector <index_l>&vl) const
{
	vp_3.clear();
	vl.clear();
	vector<edge> ve;
	vector<int> vet_same;			//用来记录每个边出现了几次与ve[i]同向的情况
	vector<int> vet_revise;			//用来记录每个边出现了几次与ve[i]逆向的情况
	for (int i(0);i!=vt.size ();++i)
	{
		int n_on(0), n_pos(0), n_neg(0);
		switch (pl.oriented_side(vp[vt[i].n[0]]))
		{
		case CGAL::ON_ORIENTED_BOUNDARY: n_on++; break;
		case CGAL::ON_POSITIVE_SIDE: n_pos++; break;
		case CGAL::ON_NEGATIVE_SIDE: n_neg++; break;
		}
		switch (pl.oriented_side(vp[vt[i].n[1]]))
		{
		case CGAL::ON_ORIENTED_BOUNDARY: n_on++; break;
		case CGAL::ON_POSITIVE_SIDE: n_pos++; break;
		case CGAL::ON_NEGATIVE_SIDE: n_neg++; break;
		}
		switch (pl.oriented_side(vp[vt[i].n[2]]))
		{
		case CGAL::ON_ORIENTED_BOUNDARY: n_on++; break;
		case CGAL::ON_POSITIVE_SIDE: n_pos++; break;
		case CGAL::ON_NEGATIVE_SIDE: n_neg++; break;
		} 

		if (n_pos == 0)
			continue;		//如果三角形没有在pl正侧的点，则不参与计算
		if (n_on==3 || n_pos==3)
			continue;		//如果整个三角形都在pl正侧，或在pl中，则也不参与计算

		triangle_3 tri(vp[vt[i].n[0]], vp[vt[i].n[1]], vp[vt[i].n[2]]);
		CGAL::cpp11::result_of <K::Intersect_3(K::Plane_3,K::Triangle_3)>::type
			result = CGAL::intersection(pl,tri);
		if (!result)
			continue;

		if (const K::Segment_3 *s = boost::get <K::Segment_3> (&*result))
		{
			int beg = add_set(vp_3,s->source());
			int end = add_set (vp_3,s->target ());
			edge e;
			vector_3 v_temp (CGAL::cross_product(vector_3(*s), pl.orthogonal_vector()));
			if (v_temp * get_normal(i) < 0)		//保证e的左侧是domain的内部
			{
				e.n[0] = end; e.n[1] = beg;
			}
			else
			{
				e.n[0] = beg; e.n[1] = end;
			}

			int ne = add_set (ve,e);
			if (vet_same.size () < ve.size ())		//如果这个edge是第一次被添加
			{
				vet_same.push_back (0);
				vet_revise.push_back (0);
			}
			if (ve[ne].n[0] == e.n[0])
				vet_same[ne]++;		//又出现了一次和ve[ne]同向的边
			else
				vet_revise[ne]++;	//又出现了一次和ve[ne]逆向的边
		}
	}
	for (int i(0);i!=ve.size ();++i)
	{
		if (vet_same[i] == 1 && vet_revise[i] == 0)		//normal situation
			continue;
		if (vet_same[i] != 1)					//因为假设stl没有非流行边界，所以每个边顶多出现两次，正反最多一次
			throw logic_error ("plane_intersect(), wrong stl format. may be nonmanifold edges");			
	}

	//如果某个边出现了正反各出现一次，那么其不参与生成loop
	int last_p(-1);		//记录目前loop最后一个edge，其实可以用一个bool类型替代。
	int count(0);
	for (int i(0);i!=ve.size ();++i)
		if (vet_revise[i] == 1)
			count++;

	while(1)
	{
		if (last_p == -1)
		{
			if (count == ve.size ())
				break;		//每个边都被添加了，皆大欢喜
			else
			{
				for (int i(0);i!=ve.size ();++i)		//寻找一个新的开始
					if (vet_revise[i] != 1)
					{
						vl.push_back(index_l());
						vl.back ().push_back (ve[i].n[0]);
						vl.back ().push_back (ve[i].n[1]);
						last_p = i;
						vet_revise[i] = 1;
						count++;
						break;
					}
			}
		}
		if (count == ve.size ())
			throw logic_error ("plane_intersection, loop is not finished");		//loop还没生成完成，边就没了

		vector<edge> proten_e;
		vector<int> e_index;
		for (int i(0);i!=ve.size ();++i)
		{
			if (vet_revise[i] == 1)
				continue;
			if (ve[i].n[0] == ve[last_p].n[1])
			{
				proten_e.push_back (ve[i]);
				e_index.push_back (i);
			}
		}
		if (proten_e.empty ())
			throw logic_error ("plane_intersection, can't find next edge");

		right_most_compare_3 temp_comp (pl,vp_3,ve[last_p]);
		int right_most_n (0);
		for (int i(0);i!=proten_e.size ();++i)
			if (temp_comp (proten_e[right_most_n],proten_e[i]))
				right_most_n = i;
		
		last_p = e_index[right_most_n];
		vet_revise[last_p] = 1;
		count++;

		if (ve[last_p].n[1] == vl.back ()[0])
			last_p = -1;
		else
			vl.back ().push_back (ve[last_p].n[1]);
	}

	//vector<index_l> vl_temp;
	//vl_temp.reserve (vl.size ()+2);
	//for (int i(0);i!=vl.size ();++i)
	//	split_index_l(vl[i],vl_temp);
	//vl.swap (vl_temp);		//将每个loop都变成简单的loop
	return 0;

}

//int stl_model::plane_intersect (const plane_3 & pl, polygon_set_2 &ps) const
//{
//	vector<point_2> vp;
//	vector<loop> vl;
//	if (plane_intersect(pl,vp,vl) != 0)
//		return 1;
//	ps.clear ();
//	
//	vector<polygon_2> vpoly;
//	vector<int> poly_to_loop;		//记录每个poly对应哪个loop
//
//	for (int i(0);i!=vl.size ();++i)
//	{
//		vpoly.push_back (polygon_2 ());
//		poly_to_loop.push_back (i);
//		vl[i].to_polygon_2(vp,vpoly.back ());
//		if (vpoly.back ().size () <=2)
//		{
//			vpoly.pop_back();		//忽略退化的边
//			poly_to_loop.pop_back();
//		}
//	}
//
//	vector<polygon_with_holes_2> vpwh;
//	vpwh.reserve (vpoly.size ());
//	vector<int> poly_ori (vpoly.size ());	//逆时针为1顺时针为2
//	for (int i(0);i!=vpoly.size ();++i)
//	{
//		if (vpoly[i].is_counterclockwise_oriented())
//			poly_ori[i] = 1;
//		else
//		{
//			poly_ori[i] = 2;
//			vpoly[i].reverse_orientation();
//		}
//	}
//
//	vector<vector<int>> vcontained (vpoly.size ());		//记录每个多边形包含哪些多边形
//	for (int i(0);i!=vpoly.size ();++i)
//	{
//		for (int j(0);j!=vpoly.size ();++j)		//查看i应该包含在那个外表面中
//		{
//			if (i==j) 
//				continue;
//			if (polygon_in_polygon_2(vp,vl[poly_to_loop[j]], vl[poly_to_loop[i]]))
//				vcontained[j].push_back (i);
//			continue;
//		}
//	}
//
//	vector<vector<int>> res;				//记录最终结果
//	res.reserve (vpoly.size ());
//	vector<int> poly_to_res (vpoly.size());		//记录每个多边形被记录在哪个最终结果中
//	for (int i(0);i!=vpoly.size ();++i)
//	{
//		if (poly_ori[i] ==2)
//			continue;
//		res.push_back(vector<int>());
//		res.back ().push_back (i);
//		poly_to_res[i] = (int)res.size ()-1;
//		res.back ().insert (res.back ().end (), vcontained[i].begin (),vcontained[i].end ());//将逆时针旋转的多边形包含的多边形都添加进去
//	}
//	for (int i(0); i!=res.size ();++i)
//	{
//		for (int j(1);j<res[i].size ();++j)
//		{
//			if (poly_ori[res[i][j]] == 1)
//			{
//				for (int k(0);k!=res[poly_to_res[res[i][j]]].size ();++k)
//				{
//					res[i].erase (remove(res[i].begin (), res[i].end (),res[poly_to_res[res[i][j]]][k]));
//				}
//			}
//
//		}
//	}
//	return 0;
//}

void split_loop (const loop& l, vector<loop> &vl)
{
	if (is_simple (l.vpi))
		vl.push_back(l);
	else
	{
		int p_begin(0),p_end(0);
		int i(0);
		for (;i!=l.vpi.size ();++i)
		{
			int i1(i+1);
			for (;i1<(int)l.vpi.size ();++i1)
				if (l.vpi[i] == l.vpi[i1])
				{
					p_begin = i;
					p_end =i1;
					break;
				}
			if (i1 != l.vpi.size ())
				break;
		}
		if (i == l.vpi.size ())
			throw logic_error ("split_loop, p is simple");
		loop l1,l2;
		int e_prev1 (p_begin), e_prev2 (p_end);
		for (int i(p_begin); i!=p_end;e_prev1 = i++)
		{
			l1.vpi.push_back (l.vpi[i]);
			if (i == p_begin)
				continue;
			l1.vei.push_back (l.vei [e_prev1]);
			l1.ve_ori .push_back (l.ve_ori[e_prev1]);
		}
		l1.vei.push_back (l.vei[e_prev1]);
		l1.ve_ori .push_back (l.ve_ori[e_prev1]);

		for (int i (p_end); i!=l.vpi.size();e_prev2 = i++)
		{
			l2.vpi.push_back (l.vpi[i]);
			if (i == p_end)
				continue;
			l2.vei.push_back (l.vei[e_prev2]);
			l2.ve_ori.push_back (l.ve_ori [e_prev2]);
		}
		l2.vei.push_back (l.vei[e_prev2]);
		l2.ve_ori.push_back (l.ve_ori [e_prev2]);

		for (int i(0);i!=p_begin;++i)
		{
			l2.vpi.push_back (l.vpi[i]);
			l2.vei.push_back (l.vei[i]);
			l2.ve_ori.push_back (l.ve_ori[i]);
		}

		split_loop (l1,vl);
		split_loop (l2,vl);
	}
}

int split_seg (std::vector<segment_3> &vseg, std::vector<bbox_3> &vbox_seg, const point_3 & po)
{
	size_t vseg_size (vseg.size ());
	bbox_3 temp_bbox (point2bbox(po));

	for (int i(0); i!=vseg_size;++i)
	{
		if (!CGAL::do_overlap(temp_bbox,vbox_seg[i]))
			continue;
		if (!vseg[i].has_on (po))
			continue;
		if (vseg[i].target ()== po || vseg[i].source ()==po)
			break;
		vseg.push_back (segment_3 (po,vseg[i].target ()));
		vbox_seg.push_back (seg2bbox(vseg.back ()));
		vseg[i] = segment_3(vseg[i].source (), po);
		vbox_seg[i] = seg2bbox(vseg[i]);
		break;
	}
	return 0;
}

void out_face_STL (std::ostream &outfile, const std::vector <point_3> &vpo, const face &fa, const plane_3 &pl,bool flag ,const std::string &s)
{
	if (flag)
		outfile<<"solid FA"<<s<<endl;

	typedef struct
	{
		int nesting_level;
	}face_info;

	using CGAL::to_double;
	typedef CGAL::Triangulation_vertex_base_2<K>					Vb;
	typedef CGAL::Triangulation_face_base_with_info_2<face_info,K>	Fbb;
	typedef CGAL::Constrained_triangulation_face_base_2 <K, Fbb>	Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb,Fb>				TDS;
	typedef CGAL::Exact_intersections_tag							TAG;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,TAG>	CDT;

	outfile.precision(20);
	CDT cdt;
	vector_3 normal (pl.orthogonal_vector());

	double a(to_double (normal.x ()));
	double b(to_double (normal.y ()));
	double c(to_double (normal.z ()));
	double d(sqrt(a*a+b*b+c*c));
	a/=d; b/=d; c/=d;

	for (int i(0);i!=fa.vl.size ();++i)
		insert_polygon (cdt,vpo,pl,fa.vl[i]);
	mark_domains (cdt);
	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end();++fit)
	{
		if (fit->info().nesting_level%2 !=1)
			continue;
		vector<point_3> vp;
		vp.push_back (pl.to_3d(fit->vertex(0)->point()));
		vp.push_back (pl.to_3d(fit->vertex(1)->point()));
		vp.push_back (pl.to_3d(fit->vertex(2)->point()));
		if (normal*oriented_normal(vp) < 0)
			reverse(vp.begin (),vp.end());
		outfile<<"  facet normal "<<a<<' '<<b<<' '<<c<<endl;
		outfile<<"    outer loop"<<endl;
		outfile<<"      vertex "<<vp[0].x()<<' '<<vp[0].y()<<' '<<vp[0].z()<<endl;
		outfile<<"      vertex "<<vp[1].x()<<' '<<vp[1].y()<<' '<<vp[1].z()<<endl;
		outfile<<"      vertex "<<vp[2].x()<<' '<<vp[2].y()<<' '<<vp[2].z()<<endl;
		outfile<<"    endloop"<<endl;
		outfile<<"  endfacet"<<endl;
	}
	if (flag)
		outfile<<"endsolid FA"<<s<<endl;
}

std::pair<int,int> intersection (const poly_frac &pf1, const plane_3 &pl1,
								 const poly_frac &pf2, const plane_3 &pl2,
								 vector<segment_3> & rese, vector<point_3> &resp,
								 std::vector <bbox_3> &vebox)
{
	
	CGAL::cpp11::result_of <K::Intersect_3(plane_3,plane_3)>:: type
		res = CGAL::intersection (pl1,pl2);
	size_t rese_size(rese.size ()), resp_size(resp.size ());
	if (res)
	{
		if (const line_3* l = boost::get <line_3> (&*res))
		{
			point_3 p(l->point());
			vector_3 v (l->to_vector());
			vector<FT> pf_seg;		//记录两个pf和这个直线的交点
			pf_seg.reserve (pf1.vl.size () * pf2.vl.size ());
			vector<bool> pf1_seg_flag, pf2_seg_flag;	//记录每两个交点之间的交线段是否在这两个pf上
			intersection_points(pf1,p,v,pf_seg);
			intersection_points(pf2,p,v,pf_seg);
			if (pf_seg.empty ())
				return pair<int,int> (0,0);		//没有相交的部分

			sort(pf_seg.begin (), pf_seg.end ());
			pf_seg.resize (unique(pf_seg.begin (),pf_seg.end()) - pf_seg.begin ());
			mark_segments(pf1,pl1,p,v,pf_seg,pf1_seg_flag);
			mark_segments(pf2,pl2,p,v,pf_seg,pf2_seg_flag);

			if (pf_seg.size () == 1)		//只有一个交点
			{
				if (point_in_polygon_3(pf1.vl,pl1,p+pf_seg[0]*v) !=2 && point_in_polygon_3 (pf2.vl,pl2,p+pf_seg[0]*v) !=2)
					resp.push_back (p+pf_seg[0]*v);
			}
			else
			{
				for (int i(0);i!=pf_seg.size ();++i)
				{
					int i1(i+1);
					if (i1 >= pf_seg.size ())
						break;
					if (pf1_seg_flag[i] && pf2_seg_flag[i])	//如果这个线段在pf1上也在pf2上
						//rese.push_back ( segment_3(p+pf_seg[i]*v, p+pf_seg[i1]*v));
						add_set_seg_3 (rese, segment_3(p+pf_seg[i]*v, p+pf_seg[i1]*v),vebox);
				}

				for (int i(0);i!=pf_seg.size ();++i)
				{
					bool pf1_back, pf1_front, pf2_back, pf2_front;		//记录这个点在pf1和pf2上是否属于某个线段
					if (i == 0)
					{ pf1_back = pf2_back = false; pf1_front = pf1_seg_flag[i]; pf2_front = pf2_seg_flag[i];}
					else if (i == pf_seg.size ()-1)
					{ pf1_back = pf1_seg_flag[i-1]; pf2_back = pf2_seg_flag[i-1]; pf1_front = pf2_front = false;}
					else
					{ 
						pf1_back = pf1_seg_flag[i-1]; pf2_back = pf2_seg_flag[i-1];
						pf1_front = pf1_seg_flag[i]; pf2_front = pf2_seg_flag[i];
					}
					if ( (!pf1_back) && (!pf1_front) && (!pf2_back) && (!pf2_front)
						&& point_in_polygon_3(pf1.vl,pl1,p+pf_seg[0]*v) !=2
						&& point_in_polygon_3(pf2.vl,pl2,p+pf_seg[0]*v) !=2)		//这个点不在任何一个线段上，并且在pf1和pf2上
						resp.push_back (p+pf_seg[i]*v);
				}
			}	
		}
	}
	return pair<int,int> ((int) (rese.size ()-rese_size), (int) (resp.size ()-resp_size));
}

void intersection_points (const std::vector<point_3> & loop, const point_3&p, const vector_3& v, std::vector<FT> &res)
{
	vector<FT> vt;
	vt.reserve (loop.size ());
	line_3 l(p,v);
	for (int i(0);i!=loop.size ();++i)
	{
		int i1(i+1);
		if (i1>=loop.size ())
			i1 = 0;
		CGAL::cpp11::result_of <K::Intersect_3(line_3,segment_3)>::type
			intersect = CGAL::intersection(l,segment_3(loop[i],loop[i1]));
		if (intersect)
		{
			if (const point_3* inter_p = boost::get<point_3> (&*intersect))	//交点是一个点
				vt.push_back (get_para (p,v,*inter_p));
			else		//交点是一个线段
			{
				vt.push_back (get_para (p,v,loop[i]));
				vt.push_back (get_para(p,v,loop[i1]));
			}
		}
	}

	sort(vt.begin (),vt.end ());
	unique_copy(vt.begin (),vt.end(), back_inserter(res));
}

void intersection_points (const std::vector<std::vector<point_3>> &vl, const point_3&p, const vector_3& v, std::vector<FT> &res)
{
	vector<FT> vt;
	vt.reserve (vl.size ()*3);
	line_3 l(p,v);
	for (int li(0);li!=vl.size ();++li)
	{
		const vector<point_3> &loop (vl[li]);
		for (int i(0);i!=vl[li].size ();++i)
		{
			int i1(i+1);
			if (i1>=loop.size ())
				i1 = 0;
			CGAL::cpp11::result_of <K::Intersect_3(line_3,segment_3)>::type
				intersect = CGAL::intersection(l,segment_3(loop[i],loop[i1]));
			if (intersect)
			{
				if (const point_3* inter_p = boost::get<point_3> (&*intersect))	//交点是一个点
					vt.push_back (get_para (p,v,*inter_p));
				else		//交点是一个线段
				{
					vt.push_back (get_para (p,v,loop[i]));
					vt.push_back (get_para(p,v,loop[i1]));
				}
			}
		}
	}
	sort(vt.begin (),vt.end ());
	std::unique_copy(vt.begin (),vt.end(), back_inserter(res));
}

void mark_segments (const poly_frac &pf, const plane_3 &pl, const point_3&p, const vector_3&v, const std::vector <FT> &segs, std::vector <bool> &res)
{
	res.resize (segs.size ());
	for (int i(0);i!=segs.size ();++i)
	{
		int i1(i+1);
		if (i1 >=segs.size ())
			break;
		if (point_in_polygon_3(pf.vl,pl,(p+((segs[i]+segs[i1])/2)*v)) ==2)
			res[i] = false;
		else
			res[i] = true;
	}
}

void disc2pf (const point_3& center_, const FT &r, const plane_3 &pl, poly_frac & res, const int& n)
{
	point_3 center (pl.projection (center_));
	res.clear();
	res.vl.push_back (poly_frac::f_loop());
	res.vl[0].reserve (n);
	vector_3 xv (pl.base1 ()), yv (pl.base2 ());
	xv = (1/CGAL::sqrt (CGAL::to_double (xv.squared_length())))*xv;
	yv = (1/CGAL::sqrt (CGAL::to_double (yv.squared_length())))*yv;
	
	double a(0);
	double da (M_PI_2/n);
	for (int i(0);i!=n;++i)
	{
		res.vl[0].push_back (center + cos(a)*r*xv + sin(a)*r*yv);
		a+=da;
	}
}


}