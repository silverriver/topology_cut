#include "geometry.h"
#include "block_system.h"
#include <iostream>
#include <fstream>
#include <CGAL\Timer.h>
using namespace std;
int main ()
{
	TC::block_system bs;
	ifstream infile ("domain.stl");
	if (!infile)
	{
		cout<<"can't find domain.stl"<<endl;
		return 0;
	}
	else
		cout<<"open domain.stl success"<<endl;

	if (bs.init_domain(infile) !=0)
	{
		cout<<"Init. block system fail"<<endl;
		return 0;
	}
	else
		cout<<"Init. block system success"<<endl;

	infile.close();
	infile.open ("fracture.dat");

	if (!infile)
	{
		cout<<"can't find fracture.dat"<<endl;
		return 0;
	}
	else
		cout<<"Open fracture.dat success"<<endl;

	if (bs.disc_frac_parser (infile,0,30) !=0)
	{
		cout<<"Parser fractures fail"<<endl;
		return 0;
	}
	else
		cout<<"Parser frastures success"<<endl;

	CGAL::Timer timer;
	timer.reset ();
	timer.start();
	cout<<bs.identify ()<<endl;
	timer.stop();
	cout<<"Time:"<< timer.time()<<endl;
	cout<<"vseg.size:"<<bs.ved .size ()<<endl;
	cout<<"bs.vfa.size ()"<<bs.vfa.size ()<<endl;
	cout<<"bs.vblk.siz ()"<<bs.vblk.size ()<<endl;
	ofstream outfile ("faces.stl");
	for (int i(0);i!=bs.vfa.size ();++i)
		bs.output_fac(outfile,i);
	outfile.close ();
	outfile.open ("blocks.stl");
	for (int i(0);i!=bs.vblk.size ();++i)
		bs.output_blk(outfile,i);
	outfile.close();

	outfile.open ("outer_bound.stl");
	for (int i(0);i!=bs.vblk_ext.size ();++i)
		bs.output_outer_bounds(outfile, i);
	outfile.close ();

	cout<<"bs.check_face_loop_edge()"<<bs.check_face_loop_edge()<<endl;
	outfile.open ("line.txt");
	bs.output_ved(outfile);
	outfile.close();

	TC::FT tot_vol (0);
	for (int i(0);i!=bs.vblk.size ();++i)
		tot_vol += bs.cal_vol(bs.vblk[i]);
	cout<<"tot_vol"<<tot_vol<<endl;



	return 0;
}