#!/usr/bin/perl
use List::Util qw(sum);
use Statistics::Descriptive;
use Math::Trig;
use Carp;
sub swrite {
    croak "usage: swrite PICTURE ARGS" unless @_;
    my $format = shift;
    $^A = "";
    formline($format,@_);
    return $^A;
}

#Argument verification starts here

$arg=join(',',@ARGV);
if($arg=~m/-h/)
	{
	print "ASSP: version 1.0\n\narguments:\n\n-i\t[*.pdb] Input file, must be of pdb file type. PDB file format for ATOM record is: (A11,a3,1x,2(a3,1x,a1,1x,i4,2x),33x,4X)\n";
	print "-fig\t[Y/N] Pictorial view of the secondary structure assignment should be generated or not (file format is PNG). default is YES\n";
	print "-h\tHelp for ASSP\n";
	print "-v\tVersion of the program in use\n";
	print "-auth\tPrasun Kumar & Manju Bansal\n";
	print "-c\tFor any querry please contact Prof. Manju Bansal (mb@mbu.iisc.ernet.in)\n";
	exit(0);
	}
if($arg=~m/-v/)
	{
	print "Version 1.0\n";
	exit(0);
	}
if($arg=~m/-c/)
	{
	print "This program is published in, if you are using this program, please cite:\n";
	exit(0);
	}
for($i=0;$i<12;$i++)
	{
	if($ARGV[$i] eq '-i')
		{
		if($ARGV[$i+1]=~m/\.pdb/)
			{

			$file=$ARGV[$i+1];
			print "Input file is $file\n\n";
			@split=split(/\./,$file);
			}
		else
			{
			print "File is not having PDB extension\nAbrupt termination.....\n\n";
			exit(0);
			}
		}
	if($ARGV[$i] eq '-fig')
		{
		$fig=lc($ARGV[$i+1]);

		if(($fig eq 'y')||($fig eq 'yes'))
			{
			print "Pictorial view of the secondary structure assignment will be generated and saved as $split[0].png\n\n";
			$red=1;
			}
		elsif(($fig eq 'n')||($fig eq 'no'))
			{
			$red=2;
			}
		else
			{
			print "Pictorial view of the secondary structure assignment will be generated and saved as $split[0].png\n\n";
			$red=2;
			}
		}

	}

#Argument verification is over
%tcode=('ALA','A','CYS','C','ASP','D','GLU','E','PHE','F','GLY','G','HIS','H','ILE','I','LYS','K','LEU','L','MET','M','ASN','N','PRO','P','GLN','Q','ARG','R','SER','S','THR','T','VAL','V','TRP','W','TYR','Y');
%ocode=('A','ALA','C','CYS','D','ASP','E','GLU','F','PHE','G','GLY','H','HIS','I','ILE','K','LYS','L','LEU','M','MET','N','ASN','P','PRO','Q','GLN','R','ARG','S','SER','T','THR','V','VAL','W','TRP','Y','TYR');
%full=('A','AlphaHelix_1','I','PiHelix_3','G','ThreeHelix_5','a','AlphaHelix_6','i','PiHelix_11','g','ThreeHelix_12','P','PolyPro_10');
$protein=$split[0];
$result=$split[0].'_nh.out';# file will be having all the parameters used for assignment
$diff_vth=$split[0].'_diff.out';# file will be having differences of the two consecutive steps
$ss_assign=$split[0].'_cont.out';#file will be having continuous stretches
$assignment=$split[0].'_assp.out';#file with assignment done by ASSIGN
#opening the PDB file for ATOM/HETATM record. For modified residues there  must be MODRES record in the PDB file
open(A,$file);#opening the PDB file
while($seq=<A>)
	{
	if(substr($seq,0,6) eq 'MODRES')
		{
		$mod[$mo]=$seq;
		$mo++;
		}
	if(((substr($seq,0,4) eq 'ATOM')||(substr($seq,0,6) eq 'HETATM'))&&(substr($seq,13,2) eq 'CA')&&($a ne substr($seq,21,5)))
		{
		$a=substr($seq,21,5);
		if(substr($seq,21,1) eq ' ')
			{
			substr($seq,21,1)='X';
			}
		if(substr($seq,0,6) eq 'HETATM')
			{
			for($md=0;$md<$mo;$md++)
				{
				if((substr($mod[$md],16,1) eq substr($seq,21,1))&&(substr($mod[$md],18,4) eq substr($seq,22,4)))
					{
					substr($seq,17,3) = substr($mod[$md],24,3);
					last;
					}
				}
			}
		$residue=substr($seq,22,4)/1;$chain=substr($seq,21,1);
		$res12=substr($seq,17,3);
		$array[$i11]=$seq;
		$reschain=$residue.$chain;
		$aminoacid{$reschain}=$res12;
		$l1[$i11]=$tcode{$res12};
		$i11++;
		}
	if(substr($seq,0,3) eq 'END')
		{
		last;
		}
	}
close(A);
splice(@mod,0,$mo);#splicing the MODRES record
$mo=0;

$l=0;
$length=@array;
open(D,">$result");#opening the *_tab.out file for writing the local parameters
open(E,">$diff_vth");#opening the *_diff.out file for writting the differences between the parameters of two consecutive steps
open(F,">$ss_assign");#opening the *_cont.out file for writing continuous stretches
open(C,">$assignment");#openning file for writing assignment done by ASSIGN
$count=0;
#parameters calculation starts from here
for($i=0;$i<$length-3;$i++)
	{
	$count++;
	#X coordinates of CA i-i+4
	$x1=substr($array[$i],30,8);
	$x2=substr($array[$i+1],30,8);
	$x3=substr($array[$i+2],30,8);
	$x4=substr($array[$i+3],30,8);
	#Y coordinates of CA i-i+4
	$y1=substr($array[$i],38,8);
	$y2=substr($array[$i+1],38,8);
	$y3=substr($array[$i+2],38,8);
	$y4=substr($array[$i+3],38,8);
	#Z coordinates of CA i-i+4
	$z1=substr($array[$i],46,8);
	$z2=substr($array[$i+1],46,8);
	$z3=substr($array[$i+2],46,8);
	$z4=substr($array[$i+3],46,8);

	$res1=substr($array[$i],22,4);
	$res2=substr($array[$i+1],22,4);
	$res3=substr($array[$i+2],22,4);
	$res4=substr($array[$i+3],22,4);

	$chain_id=substr($array[$i],21,1);#chain id
	$x12=$x2-$x1;
	$y12=$y2-$y1;
	$z12=$z2-$z1;

	$x23=$x3-$x2;
	$y23=$y3-$y2;
	$z23=$z3-$z2;

	$x34=$x4-$x3;
	$y34=$y4-$y3;
	$z34=$z4-$z3;

	$x13=$x3-$x1;
	$y13=$y3-$y1;
	$z13=$z3-$z1;

	$x24=$x4-$x2;
	$y24=$y4-$y2;
	$z24=$z4-$z2;

	$dist12=sqrt($x12**2+$y12**2+$z12**2);
	$dist23=sqrt($x23**2+$y23**2+$z23**2);
	$dist34=sqrt($x34**2+$y34**2+$z34**2);
	$dist12=sprintf("%.1f",$dist12);
	$dist23=sprintf("%.1f",$dist23);
	$dist34=sprintf("%.1f",$dist34);

	if(($dist12<=4)&&($dist23<=4)&&($dist34<=4)&&(($res2-$res1)==1)&&(($res3-$res2)==1)&&(($res4-$res3)==1))
		{
		$x12=sprintf("%.3f",$x12);
		$y12=sprintf("%.3f",$y12);
		$z12=sprintf("%.3f",$z12);
		$x23=sprintf("%.3f",$x23);
		$y23=sprintf("%.3f",$y23);
		$z23=sprintf("%.3f",$z23);
		$x34=sprintf("%.3f",$x34);
		$y34=sprintf("%.3f",$y34);
		$z34=sprintf("%.3f",$z34);
		$divx13=$x12-$x23;
		$divy13=$y12-$y23;
		$divz13=$z12-$z23;

		$divx24=$x23-$x34;
		$divy24=$y23-$y34;
		$divz24=$z23-$z34;
		$divx13=sprintf("%.3f",$divx13);
		$divy13=sprintf("%.3f",$divy13);
		$divz13=sprintf("%.3f",$divz13);
		$divx24=sprintf("%.3f",$divx24);
		$divy24=sprintf("%.3f",$divy24);
		$divz24=sprintf("%.3f",$divz24);

		$div13=sqrt($divx13**2+$divy13**2+$divz13**2);
		$div24=sqrt($divx24**2+$divy24**2+$divz24**2);

		$c1=$divy13*$divz24-$divz13*$divy24;
		$c2=$divz13*$divx24-$divx13*$divz24;
		$c3=$divx13*$divy24-$divy13*$divx24;
		$c11[$i]=$c1/sqrt($c1*$c1+$c2*$c2+$c3*$c3);
		$c22[$i]=$c2/sqrt($c1*$c1+$c2*$c2+$c3*$c3);
		$c33[$i]=$c3/sqrt($c1*$c1+$c2*$c2+$c3*$c3);

		$ang2=57.2958*acos(($x24*$x13+$y24*$y13+$z24*$z13)/(sqrt($x24**2+$y24**2+$z24**2)*sqrt($x13**2+$y13**2+$z13**2)));
		$ang2=sprintf("%.3f",$ang2);

		$nturn=1;
		$height=($x23*$c1+$y23*$c2+$z23*$c3)/(sqrt($c1*$c1+$c2*$c2+$c3*$c3));
		if($height<0.0001)
			{
			$nturn=-1;
			$height=$height*(-1);#height calculation
			}
		$twist=($divx13*$divx24+$divy13*$divy24+$divz13*$divz24)/(sqrt($divx13*$divx13+$divy13*$divy13+$divz13*$divz13)*sqrt($divx24*$divx24+$divy24*$divy24+$divz24*$divz24));

		$costheta=1-$twist;
		$radius=sqrt($div13*$div24)/(2*$costheta);

		$twist=acos($twist);
		$twist =57.2958*$twist;#twist calculation

		$xd13=$divx13/$div13;$yd13=$divy13/$div13;$zd13=$divz13/$div13;
		$xd24=$divx24/$div24;$yd24=$divy24/$div24;$zd24=$divz24/$div24;
		$ox1=$x2-($radius*$xd13);$oy1=$y2-($radius*$yd13);$oz1=$z2-($radius*$zd13);
		$ox2=$x3-($radius*$xd24);$oy2=$y3-($radius*$yd24);$oz2=$z3-($radius*$zd24);

		$dox=$ox2-$ox1;$doz=$oz2-$oz1;$doy=$oy2-$oy1;
		$dox1=$dox/sqrt($dox*$dox+$doy*$doy+$doz*$doz);
		$doy1=$doy/sqrt($dox*$dox+$doy*$doy+$doz*$doz);
		$doz1=$doz/sqrt($dox*$dox+$doy*$doy+$doz*$doz);
		$dox1=sprintf("%.3f",$dox1);
		$doy1=sprintf("%.3f",$doy1);
		$doz1=sprintf("%.3f",$doz1);

		$height1=($x23*$dox+$y23*$doy+$z23*$doz)/(sqrt($dox**2+$doy**2+$doz**2));
		if($height1<0.0001)
			{
			$height1=$height1*(-1);#height calculation
			}

		if($nturn==-1)
			{
			$twist=360-$twist;
			}
		$angle=57.2958*acos(($dox1*$c11[$i]+$doy1*$c22[$i]+$doz1*$c33[$i])/((sqrt($c11[$i]**2+$c22[$i]**2+$c33[$i]**2))*sqrt($dox1**2+$doy1**2+$doz1**2)));
		if(($angle>35)&&($angle<145)&&($ang2<10))
			{
		#	print "prasun\t$height\t$height1\n";
			$height=$height1;
			}

		$a=(($y12*$z23)-($y23*$z12));
		$b=(($z12*$x23)-($x12*$z23));
		$c=(($x12*$y23)-($y12*$x23));

		$a1=(($y23*$z34)-($y34*$z23));
		$b1=(($z23*$x34)-($x23*$z34));
		$c1=(($x23*$y34)-($y23*$x34));
		$cos=((($a*$a1)+($b*$b1)+($c*$c1))/((sqrt(($a*$a)+($b*$b)+($c*$c)))*(sqrt(($a1*$a1)+($b1*$b1)+($c1*$c1)))));
		$theta=acos($cos);
		if( $theta > 0.0 && abs($theta - 1.0) < 0.000001 )
			{
			$theta -= 0.000001;
			}
		elsif( ($theta < 0.0) && abs($theta + 1.0) < 0.000001 )
			{
			$theta += 0.000001;
			}
		$theta =57.2958*$theta;#vtor calculation
		$stp=$a1*$x12+$b1*$y12+$c1*$z12;
		if($stp<0)
			{
			$theta=(-1*$theta)+360;
			}
		else
			{
			$theta=$theta;
			}
		$xd13=$divx13/$div13;$yd13=$divy13/$div13;$zd13=$divz13/$div13;
		$xd24=$divx24/$div24;$yd24=$divy24/$div24;$zd24=$divz24/$div24;

		$twist=sprintf("%.2f",$twist);
		$height=sprintf("%.2f",$height);
		$theta=sprintf("%.2f",$theta);
		$radius=sprintf("%.2f",$radius);
		$hi=$i1+1;
		if($string ne '')
			{
			@split_string=split(' ',$string);
			if(($hi-$split_string[0])==1)
				{
				$td=abs($split_string[10]-$twist);
				$hd=abs($split_string[11]-$height);
				$vd=abs($split_string[12]-$theta);
				$td=sprintf("%.2f",$td);
				$hd=sprintf("%.2f",$hd);
				$vd=sprintf("%.2f",$vd);
				$string_diff = swrite(<<'END',$split_string[9],$split_string[1],$split_string[2],'',substr($array[$i+3],22,4),$l1[$i+3],'',$split_string[10],'',$twist,'',$td,'',$split_string[11],'',$height,'',$hd,'',$split_string[12],'',$theta,'',$vd,'',$split_string[3],$split_string[4],'',substr($array[$i+2],22,4),$l1[$i+2]);
@<@>>>@>@<@>>>@>@<@>>>>>@|@>>>>>@|@>>>>>@|@<<<@|@<<<@|@<<<@|@>>>>>@|@>>>>>@|@>>>>>@|@>>>@>@<@>>>@>
END
				print E $string_diff;#writing in a file *_diff.out
				}
			}
		if($count>2)
			{
			$bangle=57.2958*acos($c11[$i-3]*$c11[$i]+$c22[$i-3]*$c22[$i]+$c33[$i-3]*$c33[$i]);#bending angle calculation
			$bangle=sprintf("%.2f",$bangle);
			$bangle1=57.2958*acos($c11[$i-1]*$c11[$i]+$c22[$i-1]*$c22[$i]+$c33[$i-1]*$c33[$i]);
			$bangle1=sprintf("%.2f",$bangle1);
			$string = swrite(<<'END' ,$hi,'',substr($array[$i],22,4),$l1[$i],'',substr($array[$i+1],22,4),$l1[$i+1],'',substr($array[$i+2],22,4),$l1[$i+2],'',substr($array[$i+3],22,4),$l1[$i+3],$chain_id,'',$twist,'',$height,'',$theta,'',$bangle,'',$radius);
@>>>@>@>>>@>@>@>>>@>@>@>>>@>@>@>>>@>@>@>@>>>>>>@>@>>>@>@>>>>>@>@>>>>>@>@>>>
END
			print D $string;#writting in file *_tab.out
			$str_nh=substr($array[$i],22,4).'_'.$chain_id;$str_nh=~s/ //g;
			$nh{$str_nh}=$string;
			}
		else
			{
			$bangle1=57.2958*acos($c11[$i-1]*$c11[$i]+$c22[$i-1]*$c22[$i]+$c33[$i-1]*$c33[$i]);
			$bangle1=sprintf("%.2f",$bangle1);
			$string = swrite(<<'END' ,$hi,'',substr($array[$i],22,4),$l1[$i],'',substr($array[$i+1],22,4),$l1[$i+1],'',substr($array[$i+2],22,4),$l1[$i+2],'',substr($array[$i+3],22,4),$l1[$i+3],$chain_id,'',$twist,'',$height,'',$theta,'','0.00','',$radius);
@>>>@>@>>>@>@>@>>>@>@>@>>>@>@>@>>>@>@>@>@>>>>>>@>@>>>@>@>>>>>@>@>>>>>@>@>>>
END
			print D $string;#writting in file *_tab.out
			$str_nh=substr($array[$i],22,4).'_'.$chain_id;$str_nh=~s/ //g;
			$nh{$str_nh}=$string;
			}
			$i1++;
		}
	else
		{
		$i1++;
		next;
		}
	}
close(D);
close(E);
#continuous stretch assignment start here TD<35; HD<1.1; VD<50;
open(E,$diff_vth);
@data_diff=<E>;
close(E);
push(@data_diff,'F    F F     F F  F       F      40    F    F    2.0   F      F      51      F  F    F  F');
$len_diff=@data_diff;

open(E,$result);
@data_tab=<E>;
close(E);
$len_tab=@data_tab;
for($i=0;$i<$len_diff;$i++)
	{
	@split_diff=split(' ',$data_diff[$i]);
	@split_diff1=split(' ',$ss[-1]);

	if(($split_diff[7]<=35)&&($split_diff[10]<=1.1)&&($split_diff[13]<=50)&&($s==0))
		{
		$ss[$s]=$data_diff[$i];
		$s++;
		}
	elsif(($split_diff[7]<=35)&&($split_diff[10]<=1.1)&&($split_diff[13]<=50)&&($s>0)&&($split_diff[0] eq $split_diff1[0])&&(($split_diff[1]-$split_diff1[1])==1))
		{
		$ss[$s]=$data_diff[$i];
		$s++;
		}
	elsif($s>0)
		{
		@split_line1=split(' ',$ss[0]);
		@split_line2=split(' ',$ss[-1]);
		$start=$split_line1[1];
		$end=($split_line2[3]-3);
		for($r=$start;$r<=$end;$r++)
			{
			$str=$r.'_'.$split_line1[0];
			@split_tab=split(' ',$nh{$str});
			$data_cont[$h]=$nh{$str};
			$height_ss[$h]=$split_tab[11];
			$twist_ss[$h]=$split_tab[10];
			$vtor_ss[$h]=$split_tab[12];
			$rad_ss[$h]=$split_tab[-1];
			$h++;
			}
		for($r1=0;$r1<$h;$r1++)
			{
			$sumt=($twist_ss[$r1]+$twist_ss[$r1+1])/2;
			$sumrad=($rad_ss[$r1]+$rad_ss[$r1+1])/2;
			$sumh=($height_ss[$r1]+$height_ss[$r1+1])/2;
			$sumt1=($twist_ss[$r1]+$twist_ss[$r1+1]+$twist_ss[$r1+2])/3;
			$sumh1=($height_ss[$r1]+$height_ss[$r1+1]+$height_ss[$r1+2])/3;
			$sumrad1=($rad_ss[$r1]+$rad_ss[$r1+1]+$rad_ss[$r1+2])/3;
			$twist_ss1=360-$twist_ss[$r1];
			$vtor_ss1=360-$vtor_ss[$r1];
			$sumh=sprintf("%.2f",$sumh);
			if($r1<$h-2)
				{
				$sumt_l=360-($twist_ss[$r1]+$twist_ss[$r1+1])/2;
				$sumt1_l=360-($twist_ss[$r1]+$twist_ss[$r1+1]+$twist_ss[$r1+2])/3;
				$str='';
				if((((($sumt>93.6)&&($sumt<103.5))||(($twist_ss[$r1]>93.6)&&($twist_ss[$r1]<103.5))||(($sumrad>2.2)&&($sumrad<2.4))||(($rad_ss[$r1]>2.2)&&($rad_ss[$r1]<2.4)))&&((($sumh>0.9)&&($sumh<2.1))||(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1))))&&($twist_ss[$r1]<180))
					{
				#	print "$height_ss[$r1]\t$height_ss[$r1+1]\t$sumh\n";
					$str .='A ';
					}
				elsif((((($sumt_l>93.6)&&($sumt_l<103.5))||(($twist_ss1>93.6)&&($twist_ss1<103.5))||(($sumrad>2.2)&&($sumrad<2.4))||(($rad_ss[$r1]>2.2)&&($rad_ss[$r1]<2.4)))&&((($sumh>0.9)&&($sumh<2.1))||(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1))))&&($twist_ss[$r1]>180))
					{
					$str .='a ';
					}
				elsif(((($rad_ss[$r1]>1.1)&&($rad_ss[$r1]<1.9))&&(($twist_ss[$r1]>223.6)&&($twist_ss[$r1]<253)))&&(($height_ss[$r1]>2.8)&&($height_ss[$r1]<=3.2)))
					{
					$str .='P ';
					}
				else
					{
					$str .='U ';
					}
				if((((($twist_ss[$r1]>103.4)&&($twist_ss[$r1]<114.9))||(($rad_ss[$r1]>2)&&($rad_ss[$r1]<=2.2)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]<180))
					{
					$str .='G ';
					}
				elsif((((($twist_ss1>103.4)&&($twist_ss1<114.9))||(($rad_ss1>2)&&($rad_ss1<=2.2)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]>180))
					{
					$str .='g ';
					}
				else
					{
					$str .='U ';
					}
				if((((($sumt1>77.9)&&($sumt1<93.6))||(($sumrad1>2.4)&&($sumrad1<=2.7))||(($twist_ss[$r1]>77.9)&&($twist_ss[$r1]<93.6))||(($rad_ss[$r1]>2.4)&&($rad_ss[$r1]<=2.7)))&&((($sumh1>0.9)&&($sumh1<2.1))||(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1))))&&($twist_ss[$r1]<180))
					{
					$str .='I ';
					}
				elsif((((($sumt1_l>77.9)&&($sumt1_l<93.6))||(($sumrad1>2.4)&&($sumrad1<=2.7))||(($twist_ss1>77.9)&&($twist_ss1<93.6))||(($rad_ss[$r1]>2.4)&&($rad_ss[$r1]<=2.7)))&&((($sumh1>0.9)&&($sumh1<2.1))||(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1))))&&($twist_ss[$r1]>180))
					{
					$str .='i ';
					}
				else
					{
					$str .='U ';
					}
				$str .=$data_cont[$r1];
				$assign[$r1]=$str;
				$str='';
				}
			if($r1==$h-2)
				{
				$sumt_l=360-($twist_ss[$r1]+$twist_ss[$r1+1])/2;
				$str='';
				if((((($sumt>93.6)&&($sumt<103.5))||(($twist_ss[$r1]>93.6)&&($twist_ss[$r1]<103.5))||(($sumrad>2.2)&&($sumrad<2.4))||(($rad_ss[$r1]>2.2)&&($rad_ss[$r1]<2.4)))&&((($sumh>0.9)&&($sumh<2.1))||(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1))))&&($twist_ss[$r1]<180))
					{
					$str .='A ';
					}
				elsif((((($sumt_l>93.6)&&($sumt_l<103.5))||(($twist_ss1>93.6)&&($twist_ss1<103.5))||(($sumrad>2.2)&&($sumrad<2.4))||(($rad_ss[$r1]>2.2)&&($rad_ss[$r1]<2.4)))&&((($sumh>0.9)&&($sumh<2.1))||(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1))))&&($twist_ss[$r1]>180))
					{
					$str .='a ';
					}
				elsif(((($rad_ss[$r1]>1.1)&&($rad_ss[$r1]<=1.9))&&(($twist_ss[$r1]>223.6)&&($twist_ss[$r1]<253)))&&(($height_ss[$r1]>2.8)&&($height_ss[$r1]<=3.2)))
					{
					$str .='P ';
					}
				else
					{
					$str .='U ';
					}
				if((((($twist_ss[$r1]>103.4)&&($twist_ss[$r1]<114.9))||(($rad_ss[$r1]>2)&&($rad_ss[$r1]<=2.2)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]<180))
					{
					$str .='G ';
					}
				elsif((((($twist_ss1>103.4)&&($twist_ss1<114.9))||(($rad_ss[$r1]>2)&&($rad_ss[$r1]<=2.2)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]>180))
					{
					$str .='g ';
					}
				else
					{
					$str .='U ';
					}
				if((((($twist_ss[$r1]>77.9)&&($twist_ss[$r1]<93.6))||(($rad_ss[$r1]>2.4)&&($rad_ss[$r1]<=2.7)))&&((($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1))))&&($twist_ss[$r1]<180))
					{
					$str .='I ';
					}
				elsif((((($twist_ss1>77.9)&&($twist_ss1<93.6))||(($rad_ss1>2.4)&&($rad_ss1<=2.7)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]>180))
					{
					$str .='i ';
					}
				else
					{
					$str .='U ';
					}
				$str .=$data_cont[$r1];
				$assign[$r1]=$str;
				$str='';
				}
			if($r1==$h-1)
				{
				$str='';
				if((((($twist_ss[$r1]>93.6)&&($twist_ss[$r1]<103.4))||(($rad_ss[$r1]>2.2)&&($rad_ss[$r1]<2.4)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]<180))
					{
					$str .='A ';
					}
				elsif((((($twist_ss1>93.6)&&($twist_ss1<103.4))||(($rad_ss[$r1]>2.2)&&($rad_ss[$r1]<2.4)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]>180))
					{
					$str .='a ';
					}
				elsif(((($rad_ss[$r1]>1.1)&&($rad_ss[$r1]<1.9))&&(($twist_ss[$r1]>223.6)&&($twist_ss[$r1]<253)))&&(($height_ss[$r1]>2.8)&&($height_ss[$r1]<=3.2)))
					{
					$str .='P ';
					}
				else
					{
					$str .='U ';
					}
				if((((($twist_ss[$r1]>103.4)&&($twist_ss[$r1]<114.9))||(($rad_ss[$r1]>2)&&($rad_ss[$r1]<=2.2)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]<180))
					{
					$str .='G ';
					}
				elsif((((($twist_ss1>103.4)&&($twist_ss1<114.9))||(($rad_ss[$r1]>2)&&($rad_ss[$r1]<=2.2)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]>180))
					{
					$str .='g ';
					}
				else
					{
					$str .='U ';
					}
				if((((($twist_ss[$r1]>77.9)&&($twist_ss[$r1]<93.6))||(($rad_ss[$r1]>2.4)&&($rad_ss[$r1]<=2.7)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]<180))
					{
					$str .='I ';
					}
				elsif((((($twist_ss1>77.9)&&($twist_ss1<93.6))||(($rad_ss1>2.4)&&($rad_ss1<=2.7)))&&(($height_ss[$r1]>0.9)&&($height_ss[$r1]<2.1)))&&($twist_ss[$r1]>180))
					{
					$str .='i ';
					}
				else
					{
					$str .='U ';
					}
				$str .=$data_cont[$r1];
				$assign[$r1]=$str;
				$str='';
				}
			}
		print F @assign;
		#print @assign;
		#print "\n";
		$stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@height_ss);
		$avg_height = $stat->mean();

		$stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@twist_ss);
		$avg_twist = $stat->mean();
		$avg_twist=sprintf("%.2f",$avg_twist);
		$avg_height=sprintf("%.2f",$avg_height);
		$split_line1[15]=$ocode{$split_line1[15]};
		$split_line2[17]=$ocode{$split_line2[17]};
		# strech counter
		$h1=$h1+1;
		# the output line after a continuous stretch stopped
		$string_cont = swrite(<<'END' ,$protein,$h1,$split_line1[15],$split_line1[0],$split_line1[14],$split_line2[17],$split_line2[0],$split_line2[16],$avg_twist,$avg_height);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>>>>>>>>>>>>>>@>>>>>
END
		print  F $string_cont;
		#print $string_cont;
		# the array of U U U
		@assign1=@assign;
		@ss=();@assign=();@height_ss=();@twist_ss=();@rad_ss=();@data_cont=();
		$s=0;$h=0;
		if(($split_diff[7]<35)&&($split_diff[10]<1.1)&&($split_diff[13]<50))
			{
			$ss[$s]=$data_diff[$i];
			$s++;
			}
		assign($string_cont,@assign1);
		}
	}
close(F);

open(F,$ss_assign);
@data_ss=<F>;
close(F);
$len_ss=@data_ss;

#function assign

sub assign
	{
	$pp2=0;
	# @_ are the passed arguments - after-stretch-output-line and array of all assigned characteristics U U U
	local($string,@pass)=@_;
	# number of groups in this stretch
	$len_pass=@pass;
	# iterate over all characteristics resp residues in this stretch - counter $m
	for($m=0;$m<$len_pass;$m++)
		{
		# first 'U'
		$alpha[$p]=substr($pass[$m],0,1);
		# second 'U'
		$three[$p]=substr($pass[$m],2,1);
		# third 'U'
		$pi[$p]=substr($pass[$m],4,1);
		# full line
		$step[$p]=$pass[$m];
		# push all info to dedicated arrays
		$p++;
		}
	# scalar yields length of array - iterate over all info for current stretch backwards
	# => removes all terminal lines which are exclusively U U U
	for($ch=(scalar@step-1);$ch>=0;$ch--)
		{
		if(($alpha[$ch] eq 'U')&&($three[$ch] eq 'U')&&($pi[$ch] eq 'U'))
			{
			splice(@step,$ch,1);
			splice(@alpha,$ch,1);
			splice(@pi,$ch,1);
			splice(@three,$ch,1);
			$ch=scalar@step;
			}
		else
			{
			last;
			}
		}
	# remove all leading lines of U U U
	for($ch=0;$ch<scalar@step;$ch++)
		{
		if(($alpha[$ch] eq 'U')&&($three[$ch] eq 'U')&&($pi[$ch] eq 'U'))
			{
			splice(@step,$ch,1);
			splice(@alpha,$ch,1);
			splice(@pi,$ch,1);
			splice(@three,$ch,1);
			$ch=-1;
			}
		else
			{
			last;
			}
		}
	$v=0;
	$lens=@step;
	# nothing to do, when stretch has length 0 or 1
	if($lens<=1)
		{
		@alpha=();@three=();@pi=();@step=();$p=0;
		next;
		}
	else
		{
		# iterate over all interesting lines of stretches
		for($n=0;$n<$lens;$n++)
			{
			# split line
			# e.g. U U I   34    34 K    35 S    36 H    37 P A    80.77  0.92   23.52   10.23  2.82
			@split=split(' ',$step[$n]);
			# i => case-insenstive matching
			if(($alpha[$n] =~m/A/i)&&($three[$n] =~m/U/)&&($pi[$n] =~m/U/))
				{
				$final[$n]='A';
				}
			elsif(($alpha[$n] =~m/P/)&&($three[$n] =~m/U/)&&($pi[$n] =~m/U/))
				{
				$final[$n]='P';
				# ?
				$pp2=1;
				}
			elsif(($alpha[$n] =~m/A/i)&&($three[$n] =~m/G/i)&&($pi[$n] =~m/U/))
				{
				# last element
				if($split[-1]<2.2)
					{
					$final[$n]='G';
					}
				else
					{
					$final[$n]='A';
					}
				}
			elsif(($alpha[$n] =~m/A/i)&&($three[$n] =~m/U/)&&($pi[$n] =~m/I/i))
				{
				if($split[-1]>2.4)
					{
					$final[$n]='I';
					}
				else
					{
					$final[$n]='A';
					}
				}
			elsif(($alpha[$n] =~m/U/)&&($three[$n] =~m/U/)&&($pi[$n] =~m/I/i))
				{
				$final[$n]='I';
				}
			elsif(($alpha[$n] =~m/U/)&&($three[$n] =~m/G/i)&&($pi[$n] =~m/U/))
				{
				$final[$n]='G';
				}
			elsif(($alpha[$n] =~m/A/i)&&($three[$n] =~m/G/i)&&($pi[$n] =~m/I/i))
				{
				if($split[-1]<=2.2)
					{
					$final[$n]='G';
					}
				elsif($split[-1]>2.4)
					{
					$final[$n]='I';
					}
				else
					{
					$final[$n]='A';
					}
				}
			elsif(($alpha[$n] =~m/U/)&&($three[$n] =~m/U/)&&($pi[$n] =~m/U/))
				{
				@split_st=split(' ',$step[$n]);
				$mod_t=360-$split_st[13];
				if(((($split_st[13]>104.25)&&($split_st[13]<=120))||(($mod_t>104.25)&&($split_st[13]>180)))&&($pp2!=1))
					{
					$final[$n]='G';
					}
				elsif((($split_st[13]<92.17)||(($mod_t<92.17)&&($split_st[13]>180)))&&($pp2!=1))
					{
					$final[$n]='I';
					}
				elsif(((($split_st[-1]>1.1)&&($split_st[-1]<1.9))||(($split_st[13]>220)&&($split_st[13]<260)))&&(($split_st[14]=>2.7)&&($split_st[14]<=3.2))&&(($split_st[13]>180)))
					{
					$final[$n]='P';
					}
				else
					{
					if(($pp2!=1)&&($split_st[13]<120))
						{
						$final[$n]='A';
						}
					else
						{
						$final[$n]='U';
						}
					}
				}
			}
		$finala='';
		$finala=join(' ',@final);
		# length of final string
		$lenf=@final;
		# number of residues in cut (!) stretch (i.e. without initial or terminal U U U)
		$lenstep=@step;
	#	print @step;
	#	print "$finala\n";
		if((($lenf<3)&&($finala=~m/P{2,}/))||(($finala!~m/P P/)&&($finala=~m/P/)))
			{
			$finala='';
			}
		if($lenf<2)
			{
			$finala='';
			}
		if(($lenf==2)&&($finala!~m/P/))
			{
			if(($finala eq 'A G')||($finala eq 'A A')||($finala eq 'G A'))
				{
				@split=split(' ',$step[0]);
				@split1=split(' ',$step[1]);
				$avg_twist=($split[13]+$split1[13])/2;
				$avg_rad=($split[-1]+$split1[-1])/2;
				$mod_t=360-$avg_twist;
				if(($avg_rad<=2.3)&&((($avg_twist>102)&&($avg_twist<=120))||(($mod_t>102)&&($avg_twist>180))))
					{
					$finala='G G';
					}
				else
					{
					$finala='';
					}
				}
			else
				{
				$finala='';
				}
			}
		if($lenf==3)
			{
			if(($finala eq 'A G G')||($finala eq 'G A G'))
				{
				$finala='G G G';
				}
			if(($finala eq 'A G A')||($finala eq 'G A A'))
				{
				$finala='A A A';
				}
			}
		if($finala ne '')
			{
			#print "$finala\n";
			@assignr= final_assign($finala,@step);

			push(@assign_raw,@assignr);
			@assignr=();
			}
		$finala='';
		@final=();@alpha=();@three=();@pi=();@step=();$p=0;$n=0;
		}
	}
#function final_assign getting called by assign
sub final_assign
	{
	# type - string (A U U I G G...), steps - residues/lines
	local($type,@steps)=@_;
	# length of the assignment string
	$length1=length($type);
	# number of residues/lines
	$steps=@steps;
	@rassign=();
	$ra=0;
	if((substr($type,0,7) eq 'A I I I')||(substr($type,0,7) eq 'G I I I')||(substr($type,0,7) eq 'I A I I')||(substr($type,0,7) eq 'A I I A'))
		{
		substr($type,0,7)='I I I I';
		}
	if((substr($type,0,7) eq 'A I A A')||(substr($type,0,7) eq 'A G A A')||(substr($type,0,7) eq 'I I A A'))
		{
		substr($type,0,7)='A A A A';
		}
	if((substr($type,0,9) eq 'A A I I I')||(substr($type,0,9) eq 'I A A I I'))
		{
		substr($type,0,9)='I I I I I';
		}
	if((substr($type,$length1-7,7) eq 'I I I A')||(substr($type,$length1-7,7) eq 'I I I G')||(substr($type,$length1-7,7) eq 'I I A I'))
		{
		substr($type,$length1-7,7)='I I I I';
		}
	if(substr($type,$length1-7,7) eq 'A A I I')
		{
		substr($type,$length1-7,7)='A A A A';
		}
	if((substr($type,0,7) eq 'A G G G')||(substr($type,0,7) eq 'G A G G')||(substr($type,0,7) eq 'G G A G')||(substr($type,0,7) eq 'A I G G'))
		{
		substr($type,0,7)='G G G G';
		}
	if((substr($type,0,5) eq 'A G G')||(substr($type,0,5) eq 'G A G'))
		{
		substr($type,0,5)='G G G';
		}
	if((substr($type,0,5) eq 'I A A')||(substr($type,0,5) eq 'G A A'))
		{
		substr($type,0,5)='A A A';
		}
	if(substr($type,$length1-5,5) eq 'G G A')
		{
		substr($type,$length1-5,5)='G G G';
		}
	if((substr($type,$length1-5,5) eq 'A A I')||(substr($type,$length1-5,5) eq 'A A G')||(substr($type,$length1-5,5) eq 'A G A')||(substr($type,$length1-5,5) eq 'A I A'))
		{
		substr($type,$length1-5,5)='A A A';
		}
	if((substr($type,$length1-7,7) eq 'G G G A')||(substr($type,$length1-7,7) eq 'G G A G'))
		{
		substr($type,$length1-7,7)='G G G G';
		}
	if(substr($type,0,7) eq 'G A A G')
		{
		substr($type,0,7)='A A A A';
		}
	for($l=0;$l<$length1;$l++)
		{
		if(substr($type,$l,9) eq 'G G A G G')
			{
			substr($type,$l,9) = 'G G G G G';
			}
		if((substr($type,$l,9) eq 'A A G A A')||(substr($type,$l,9) eq 'A A I A A'))
			{
			substr($type,$l,9) = 'A A A A A';
			}
		if((substr($type,$l,9) eq 'G G A A G')||(substr($type,$l,9) eq 'G A G A G'))
			{
			substr($type,$l,9) = 'G G G G G';
			}
		if((substr($type,$l,5) eq 'G A I')||(substr($type,$l,5) eq 'I A G')||(substr($type,$l,5) eq 'I G A')||(substr($type,$l,5) eq 'G I A')||(substr($type,$l,5) eq 'A I A')||(substr($type,$l,5) eq 'A G I')||(substr($type,$l,5) eq 'A I G'))
			{
			#print "prasun\n";
			#print $steps[$l];
			substr($type,$l,5) = 'A A A';
			}
		if(substr($type,$l,9) eq 'I I A I I')
			{
			substr($type,$l,9) = 'I I I I I';
			}
		if((substr($type,$l,7) eq 'I I A I')||(substr($type,$l,7) eq 'I I G I'))
			{
			substr($type,$l,7) = 'I I I I';
			}
		if((substr($type,$l,7) eq 'A A I U')||(substr($type,$l,7) eq 'A A G U'))
			{
			substr($type,$l,7) = 'A A A U';
			}
		}
	if((substr($type,$length1-9,9) eq 'A A A I I')||(substr($type,$length1-9,9) eq 'A A A G A'))
		{
		substr($type,$length1-9,9)='A A A A A';
		}
	if(substr($type,$length1-9,9) eq 'I I I A A')
		{
		substr($type,$length1-9,9)='I I I I I';
		}
	# individual chars
	@split_t=split(' ',$type);
	push(@split_t,'F');
	# number of chars
	$len_t=@split_t;
	#print "type\t$type\t$len_t\n";
	# for each char
	for($l1=0;$l1<$len_t;$l1++)
		{
		# char is equal next char
		if($split_t[$l1+1] eq $split_t[$l1])
			{
			$assignment[$as]=$steps[$l1];
		#	print "$split_t[$l1+1] \t $split_t[$l1]\n";
			$as++;
			$last_step=$steps[$l1+1];
			# type of the sse
			$char=$split_t[$l1];
			}
		else
			{
			if($as>0)
				{
				push(@assignment,$last_step);
				$start=$assignment[0];
				#print "$start\n";
				$end=$assignment[-1];
				@split_start=split(' ',$start);
				@split_end=split(' ',$end);
				#print "$char\t".substr($start,0,5)."\n";
				if((substr($start,0,5)=~m/[AGI]/)||((substr($start,0,5) eq 'U U U')&&($l1!=0)&&($split_start[13]<180)))
					{
					$type_h=$char;
					}
				elsif((substr($start,0,5)=~m/[agi]/)||((substr($start,0,5) eq 'U U U')&&($l1!=0)&&($split_start[13]>180)))
					{
					# to lower case
					$type_h=lc($char);
					}
				elsif(substr($start,0,5)=~m/P/)
					{
					$type_h=$char;
					}
				if((substr($start,0,5)=~m/[AGIPagi]/)||(($char=~m/[AGIagi]/)&&($l1!=0)))
					{
					$hel_type=$full{$type_h};$type_h='';
					$start_res=$ocode{$split_start[7]};
					$end_res=$ocode{$split_end[9]};
					@split_hel=split('_',$hel_type);
					$hellen=$split_end[8]-$split_start[6]+1;
					$string = swrite(<<'END' ,$split_hel[0],$start_res,$split_start[12],$split_start[6],$end_res,$split_end[12],$split_end[8],$split_hel[1],$hellen);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
					$rassign[$ra]=$string;
					$ra++;
					}
				splice(@assignment,0,$as+1);
				$as=0;
				}
			}
		}
	splice(@alpha_t,0,$a11);
	splice(@three_t,0,$t11);
	splice(@pi_t,0,$p11);
	$a11=0;$t11=0;$p11=0;
	return @rassign;
	}
#ASSIGN result reorganization
$len_raw=@assign_raw;
#print @assign_raw;
# iterate over all secondary structure elements
for($l=0;$l<$len_raw;$l++)
	{
	# current element
	@split_raw=split(' ',$assign_raw[$l]);
	$type=$split_raw[7];
	$sr=$split_raw[3];
	$er=$split_raw[6];
	$chainr=$split_raw[2];
	$length=$split_raw[8];
#previous
	@split_raw0=split(' ',$assign_raw[$l-1]);
	$type0=$split_raw0[7];
	$sr0=$split_raw0[3];
	$er0=$split_raw0[6];
	$chainr0=$split_raw0[2];
	$length0=$split_raw0[8];
#next
	@split_raw1=split(' ',$assign_raw[$l+1]);
	$type1=$split_raw1[7];
	$sr1=$split_raw1[3];
	$er1=$split_raw1[6];
	$chainr1=$split_raw1[2];
	$length1=$split_raw1[8];

	if((($type==6)||($type==10)||($type==11)||($type==12))&&((($type0!=6)||($type0!=10)||($type0!=11)||($type0!=12))||($chainr ne $chainr0)))
		{
		if(($er>=$sr1)&&($chainr eq $chainr1))
			{
			if((($type==6)&&($length<5))||(($type==11)&&($length<6))||(($type==12)&&($length<4))||(($type==10)&&($length<4)))
				{
				splice(@assign_raw,$l,1);
				$l--;
				}
			else
				{
				$er=$er-1;$er1=$er.$chainr;
				$length=$length-1;
				$string = swrite(<<'END' ,$split_raw[0],$split_raw[1],$chainr,$sr,$aminoacid{$er1},$chainr,$er,$type,$length);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
				$assign_raw[$l]=$string;
				}
			}
		if(($sr<=$er0)&&($chainr eq $chainr0))
			{
			if((($type==6)&&($length<5))||(($type==11)&&($length<6))||(($type==12)&&($length<4))||(($type==10)&&($length<4)))
				{
				splice(@assign_raw,$l,1);
				$l--;
				}
			else
				{
				$sr=$sr+1;
				$length=$length-1;$sr1=$sr.$chainr;$er1=$er.$chainr;
				$string = swrite(<<'END' ,$split_raw[0],$aminoacid{$sr1},$chainr,$sr,$aminoacid{$er1},$chainr,$er,$type,$length);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
				$assign_raw[$l]=$string;
				}
			}
		}
	elsif($length==2)
		{
		next;
		}
	elsif(($type == $type1)&&($er==$sr1)&&($chainr eq $chainr1))
		{
		# traverse all lines of diff
		for($r=0;$r<$len_tab;$r++)
			{
			@splitr=split(' ',$data_tab[$r]);
			if(($splitr[1]==$er)&&($splitr[9] eq $chainr))
				{
				if($splitr[13]< 60)
					{
					$er=$er1;
					$split_raw[4]=$split_raw1[4];
					$split_raw[8]=$length+$length1-1;
					$string = swrite(<<'END' ,$split_raw[0],$split_raw[1],$chainr,$sr,$split_raw[4],$chainr,$er,$type,$split_raw[8]);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
					$assign_raw[$l]=$string;
					splice(@assign_raw,$l+1,1);
					$l=$l-1;
					last;
					}
				else
					{
					reorganize($assign_raw[$l],$assign_raw[$l+1]);
					last;
					}
				}
			}
		}
	elsif(($type != $type1)&&($er==$sr1)&&($chainr eq $chainr1))
		{
		reorganize($assign_raw[$l],$assign_raw[$l+1]);
		}
	}
#reorganization fumction
sub reorganize
	{
	local($line_assign,$line_assign1)=@_;
	@split_line=split(' ',$line_assign);
	@split_line1=split(' ',$line_assign1);
	#print "$line_assign.$line_assign1";
	$htype1=$split_line[0];
	$start1=$split_line[3];
	$end1=$split_line[6];
	$length1=$split_line[8];
	$type_1=$split_line[7];

	$htype2=$split_line1[0];
	$start2=$split_line1[3];
	$end2=$split_line1[6];
	$length2=$split_line1[8];
	$type_2=$split_line1[7];

	if(((($htype2 eq 'AlphaHelix')&&($length2>4))||(($htype2 eq 'PiHelix')&&($length2 > 5))||(($htype2 eq 'ThreeHelix')&&($length2>3)))&&((($htype1 eq 'AlphaHelix')&&($length1>4))||(($htype1 eq 'PiHelix')&&($length1 > 5))||(($htype1 eq 'ThreeHelix')&&($length1>3)))&&((($type_1<6)&&($type_2<6))||(($type_1>5)&&($type_2>5))))
		{
		if((($htype1 eq 'PiHelix')||($htype1 eq 'ThreeHelix'))||($htype1 eq $htype2))
			{
			$start2=$start2+1;
			$length2=$length2-1;
			}
		if(($htype1 eq 'AlphaHelix')&&(($htype2 eq 'PiHelix')||($htype2 eq 'ThreeHelix')))
			{
			for($ps=0;$ps<$len_ss;$ps++)
				{
				# cont.out
				@splitss=split(' ',$data_ss[$ps]);
				@splitss1=split(' ',$data_ss[$ps+1]);
				#                                               chainId matching
				if(($splitss[8]==$end1)&&($splitss1[6]==$end1)&&($split_line[2] eq $splitss[12]))
					{
					if(((substr($data_ss[$ps],0,5)=~m/I/)&&(substr($data_ss[$ps+1],0,5)=~m/I/)&&($split_line1[0] eq 'PiHelix'))||((substr($data_ss[$ps],0,5)=~m/G/)&&(substr($data_ss[$ps+1],0,5)=~m/G/)&&($split_line1[0] eq 'ThreeHelix')))
						{
						$end1=$end1-1;
						$length1=$length1-1;
						}
					else
						{
						$start2=$start2+1;
						$length2=$length2-1;
						}
					}
				if(($splitss1[6]> $end1)&&($splitss1[12] eq $split_line[2]))
					{
					last;
					}
				}
			}
		$end11=$end1.$split_line[2];
		$string = swrite(<<'END' ,$htype1,$split_line[1],$split_line[2],$start1,$aminoacid{$end11},$split_line[2],$end1,$type_1,$length1);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$assign_raw[$l]=$string;
		$start21=$start2.$split_line1[2];
		$string = swrite(<<'END' ,$htype2,$aminoacid{$start21},$split_line1[2],$start2,$split_line1[4],$split_line1[2],$split_line1[6],$type_2,$length2);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$assign_raw[$l+1]=$string;
		}
	elsif(((($htype2 eq 'AlphaHelix')&&($length2>4))||(($htype2 eq 'PiHelix')&&($length2 > 5))||(($htype2 eq 'ThreeHelix')&&($length2>3)))&&((($htype1 eq 'AlphaHelix')&&($length1==4))||(($htype1 eq 'PiHelix')&&($length1 == 5))||(($htype1 eq 'ThreeHelix')&&($length1==3))))
		{
		$start2=$start2+1;
		$length2=$length2-1;
		$start21=$start2.$split_line1[2];
		$string = swrite(<<'END' ,$htype2,$aminoacid{$start21},$split_line1[2],$start2,$split_line1[4],$split_line1[2],$split_line1[6],$type_2,$length2);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$assign_raw[$l+1]=$string;
		}
	elsif(((($htype2 eq 'AlphaHelix')&&($length2==4))||(($htype2 eq 'PiHelix')&&($length2 == 5))||(($htype2 eq 'ThreeHelix')&&($length2==3)))&&((($htype1 eq 'AlphaHelix')&&($length1>4))||(($htype1 eq 'PiHelix')&&($length1 > 5))||(($htype1 eq 'ThreeHelix')&&($length1>3))))
		{
		$end1=$end1-1;
		$length1=$length1-1;
		$end11=$end1.$split_line[2];
		$string = swrite(<<'END' ,$htype1,$split_line[1],$split_line[2],$start1,$aminoacid{$end11},$split_line[2],$end1,$type_1,$length1);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$assign_raw[$l]=$string;
		}
	elsif(((($htype2 eq 'AlphaHelix')&&($length2==4))||(($htype2 eq 'PiHelix')&&($length2 == 5))||(($htype2 eq 'ThreeHelix')&&($length2==3)))&&((($htype1 eq 'AlphaHelix')&&($length1==4))||(($htype1 eq 'PiHelix')&&($length1 == 5))||(($htype1 eq 'ThreeHelix')&&($length1==3))))
		{
		if($type_1 == $type_2)
			{
			$end1=$end2;
			$split_line[4]=$split_line1[4];
			$length1=$length1+$length2-1;
			$lp=$l+1;
			splice(@assign_raw,$lp,1);
			$l=$l-1;
			$string = swrite(<<'END' ,$htype1,$split_line[1],$split_line[2],$start1,$split_line[4],$split_line[2],$end1,$type_1,$length1);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$assign_raw[$l]=$string;
			}
		elsif($length1>$length2)
			{
			$end1=$end2;
			$length1=$length2+$length-1;
			$end11=$end1.$split_line[2];
			$htype2=pi_twist_check($start1,$end1,$split_line[2],$length1);
			$string = swrite(<<'END' ,$htype1,$split_line[1],$split_line[2],$start1,$aminoacid{$end11},$split_line[2],$end1,$type_1,$length1);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$assign_raw[$l]=$string;
			$lp=$l+1;
			splice(@assign_raw,$lp,1);
			$l=$l-1;

			}
		elsif($length2>$length1)
			{
			$length2=$length2+$length-1;
			$start2=$start1;
			$start21=$start2.$split_line1[2];
			$htype2=pi_twist_check($start2,$end2,$split_line1[2],$length2);
			$string = swrite(<<'END' ,$htype2,$aminoacid{$start21},$split_line1[2],$start2,$split_line1[4],$split_line1[2],$split_line1[6],$type_2,$length2);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$assign_raw[$l]=$string;
			$lp=$l+1;
			splice(@assign_raw,$lp,1);
			$l=$l-1;
			}
		}
	elsif(((($htype2 eq 'AlphaHelix')&&($length2<4))||(($htype2 eq 'PiHelix')&&($length2 < 5))||(($htype2 eq 'ThreeHelix')&&($length2<3)))&&((($htype1 eq 'AlphaHelix')&&($length1>=4))||(($htype1 eq 'PiHelix')&&($length1 >= 5))||(($htype1 eq 'ThreeHelix')&&($length1>=3))))
		{
		$end1=$end2;
		$length1=$length2+$length-1;
		$end11=$end1.$split_line[2];
		$string = swrite(<<'END' ,$htype1,$split_line[1],$split_line[2],$start1,$aminoacid{$end11},$split_line[2],$end1,$type_1,$length1);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$assign_raw[$l]=$string;
		$lp=$l+1;
		splice(@assign_raw,$lp,1);
		$l=$l-1;
		}
	elsif(((($htype2 eq 'AlphaHelix')&&($length2>=4))||(($htype2 eq 'PiHelix')&&($length2 >= 5))||(($htype2 eq 'ThreeHelix')&&($length2>=3)))&&((($htype1 eq 'AlphaHelix')&&($length1<4))||(($htype1 eq 'PiHelix')&&($length1 < 5))||(($htype1 eq 'ThreeHelix')&&($length1<3))))
		{
		$length2=$length2+$length-1;
		$start2=$start1;
		$start21=$start2.$split_line1[2];
		$string = swrite(<<'END' ,$htype2,$aminoacid{$start21},$split_line1[2],$start2,$split_line1[4],$split_line1[2],$split_line1[6],$type_2,$length2);
@<<<<<<<<<<@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$assign_raw[$l]=$string;
		$lp=$l+1;
		splice(@assign_raw,$lp,1);
		$l=$l-1;
		}
	}
#	print "\n\n";print @assign_raw;
chomp(@assign_raw);
for($q=0;$q<$len_raw;$q++)
	{
	@split_final=split(' ',$assign_raw[$q]);
	# split_final[8] = length
	if($split_final[8]==2)
		{
		delete($assign_raw[$q]);
		next;
		}
	if($split_final[8] eq '')
		{
		delete($assign_raw[$q]);
		next;
		}
	if(($split_final[0] eq 'AlphaHelix')&&($split_final[8]>3))
		{
		$split_final[0]='AlphaHelix';
		# type == 'A'
		if($split_final[7]==1)
			{
			@splita=split(' ',$alphar[$ar1-1]);
			if(($pi2alpha==1)&&($splita[7]==$split_final[3]-1)&&($split_final[2] eq $splita[3]))
				{
				$split_final[3]=$splita[4];
				splice(@alphar,$ar1-1,1);
				$ar1--;
				$pi2alpha=0;
				$split_final[8]+=$splita[9];
				}
			$ar1=$ar1+1;
			$str=swrite(<<'END' ,$split_final[0],$ar1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],$split_final[7],$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$alphar[$ar1-1]=$str;
			}
        # type == 'a'
		elsif($split_final[7]==6)
			{
			$ah1=$ah1+1;
			$str=swrite(<<'END' ,$split_final[0],$ah1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],$split_final[7],$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$alphah[$ah1-1]=$str;
			}
		}
	if((($split_final[0] eq 'ThreeHelix')&&($split_final[8]>2)))
		{
		$split_final[0]='ThreeHelix';
		if(($split_final[7]==1)||($split_final[7]==5))
			{
			$split_final[7]=5;
			$tr1=$tr1+1;
			$str=swrite(<<'END' ,$split_final[0],$tr1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],$split_final[7],$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$threer[$tr1-1]=$str;
			}
		if(($split_final[7]==6)||($split_final[7]==12))
			{
			$split_final[7]=12;
			$th1=$th1+1;
			$str=swrite(<<'END' ,$split_final[0],$th1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],$split_final[7],$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$threeh[$th1-1]=$str;
			}
		}
	if(($split_final[0] eq 'PiHelix')&&($split_final[8]>4))
		{
		if($split_final[7]==3)
			{
			$split_final[0]=pi_twist_check($split_final[3],$split_final[6],$split_final[2],$split_final[8]);#print "$split_final[0]\n";
			if($split_final[0] eq 'AlphaHelix')
				{
				@splita=split(' ',$alphar[$ar1-1]);
				# split final is the line of PDB SSE annotation
				# splita should be the previous annotation
				if(($splita[7]==$split_final[3]-1)&&($split_final[2] eq $splita[3]))
					{
					$split_final[3]=$splita[4];
					splice(@alphar,$ar1-1,1);
					$ar1--;
					$pi2alpha=1;
					$split_final[8]+=$splita[9];
					}
				$ar1=$ar1+1;
				$str=swrite(<<'END' ,$split_final[0],$ar1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],1,$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
				$alphar[$ar1-1]=$str;
				}
			else
				{
				$pr1=$pr1+1;
				$str=swrite(<<'END' ,$split_final[0],$pr1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],$split_final[7],$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
				$pir[$pr1-1]=$str;
				}

			}
		elsif($split_final[7]==11)
			{
			$ph1=$ph1+1;
			$str=swrite(<<'END' ,$split_final[0],$ph1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],$split_final[7],$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$pih[$ph1-1]=$str;
			}
		}
	if($split_final[0] eq 'PolyPro')
		{
		$pp1=$pp1+1;
		$str=swrite(<<'END' ,$split_final[0],$pp1,$split_final[1],$split_final[2],$split_final[3],$split_final[4],$split_final[5],$split_final[6],$split_final[7],$split_final[8]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$pp[$pp1-1]=$str;
		}
	}
print C @alphar;
print C @alphah;
print C @threer;
print C @threeh;
print C @pir;
print C @pih;
print C @pp;
close(C);
open(D,"$result");
@data=<D>;
close(D);
open(D,">$result");
print D "============================================================================
Sl.No.  Res1    Res2    Res3    Res4 Ch  Twist     h    Vtor       BA Radius
============================================================================\n";
print D @data;
print D "================================================================================
 DESCRIPTION OF THE RECORD
================================================================================
 Res1      : CA atom of ith residue of the repeating unit.
 Res2      : CA atom of (i+1)th residue of the repeating unit.
 Res3      : CA atom of (i+2)th residue of the repeating unit.
 Res4      : CA atom of (i+3)th residue of the repeating unit.
 ch        : Chain identifier of the repeating unit.
 Twist     : Translation along the helical axis per repeating unit (Deg.).
 h         : Rise per Residue along the helical axis per repeating unit (Ang.).
 Vtor      : Virtual torsional angle per repeating unit (Deg.).
 BA        : Bending angle b/w local helix axes of two successive turns(Deg.).
 Radius    : Radius of least squares circle fitted per repeating unit (Ang.).
================================================================================";
close(D);
open(E,"$diff_vth");
@data=<E>;
close(E);
open(E,">$diff_vth");
print E "==================================================================================================
ch  Res1    Res5  Twist1  Twist2   ΔTwist  h1    h2    Δh    Vtor1   Vtor2   ΔVtor     Res2   Res4
==================================================================================================\n";
print E @data;
print E "==================================================================================================
 DESCRIPTION OF THE RECORD
==================================================================================================
 ch        : Chain identifier of the repeating unit 1 and 2.
 Res1      : CA atom of 1st residue of the repeating unit 1.
 Res5      : CA atom of 4th residue of the repeating unit 2.
 Twist#    : Translation along the helical axis for repeating unit # (Deg.).
 ΔTwist    : Absolute difference between Twist1 and Twist2 (Deg.).
 h#        : Rise per Residue (height) along the helical axis forr repeating unit # (Ang.).
 Δh        : Absolute difference between Height1 and Height2 (Ang.).
 Vtor#     : Virtual torsional angle for repeating unit # (Deg.).
 ΔVtor     : Absolute difference between Vtor1 and Vtor2 (Deg.).
 Res2      : CA atom of 2nd residue of the repeating unit 1.
 Res4      : CA atom of 3rd residue of the repeating unit 2.
==================================================================================================";
close(E);
open(F,"$ss_assign");
@data=<F>;
close(F);
open(F,">$ss_assign");
print F "===================================================================================
1 2 3  Sl.No.  Res1   Res2    Res3    Res4 Ch  Twist    h     Vtor       BA  Radius
===================================================================================\n";
print F @data;
print F "=================================================================================
 DESCRIPTION OF THE RECORD
=================================================================================
 1         : Twist, Rise per Residue and Radius satisfying the criteria of either
             alpha helix or PPII (can have A, a, P, U).
 2         : Twist, Rise per Residue and Radius satisfying the criteria of 310
             helix (can have G, g, U).
 3         : Twist, Rise per Residue and Radius satisfying the criteria of pi
             helix (can have I, i, U).
    A = Right-handed α helix                      a = Left-handed α helix
    G = Right-handed 310 helix                    g = Left-handed 310 helix
    I = Right-handed π helix                      i = Left-handed π helix
    P = PolyProline-II                            U = Undefined
 Res1      : CA atom of ith residue of the repeating unit.
 Res2      : CA atom of (i+1)th residue of the repeating unit.
 Res3      : CA atom of (i+2)th residue of the repeating unit.
 Res4      : CA atom of (i+3)th residue of the repeating unit.
 ch        : Chain identifier of the repeating unit.
 Twist     : Translation along the helical axis per repeating unit (Deg.).
 h         : Rise per Residue along the helical axis per repeating unit (Ang.).
 Vtor      : Virtual torsional angle per repeating unit (Deg.).
 BA        : Bending angle b/w local helix axes of two successive turns(Deg.).
 Radius    : Radius of least squares circle fitted per repeating unit (Ang.).
================================================================================";
close(F);

sub pi_twist_check
	{
	$tw=0;$rad=0;
	local($start,$end,$chain_p,$length_p)=@_;
	#print "$start\t$end\t$chain\t$length_p\t$aiwe\n";
	open(F,$result);
	@rtab=<F>;
	close(F);
	for($qw=0;$qw<scalar@rtab;$qw++)
		{
		@splitr=split(' ',$rtab[$qw]);
		if(($splitr[1]>=$start-1)&&($splitr[7]<=$end+1)&&($splitr[9] eq $chain_p))
			{
			#print $rtab[$qw];
			$tw+=$splitr[10]/($length_p-1);
			$rad+=$splitr[-1]/($length_p-1);
			}
		if(($splitr[7]>$end+1)&&($splitr[9] eq $chain_p))
			{
			last;
			}
		}
	#print $tw."\t$rad\n";
	if(($tw<93.8)&&($rad>2.45))
		{
		$htype='PiHelix';
		}
	elsif(($tw>102)&&($rad<2.3))
		{
		$htype='ThreeHelix';
		}
	else
		{
		$htype='AlphaHelix';
		}
	$tw=0;
	#print "$htype\n";
	return $htype;
	}
#adding for strands

$twist=199.92;$stdt=18.37;
$height=3.01;$stdh=0.44;
$vtor=201.91;$stdv=21.28;
$rad=0.95;$stdr=0.18;

$t2=$twist+2*$stdt;$t1=$twist-2*$stdt;
$h2=$height+2*$stdh;$h1=$height-2*$stdh;
$v2=$vtor+2*$stdv;$v1=$vtor-2*$stdv;
$r2=$rad+2*$stdr;$r1=$rad-2*$stdr;

open(A,$assignment);
@helix=grep {m/Helix|PolyPro/} <A>;
close(A);
for($i=0;$i<scalar@helix;$i++)
	{
	@split=split(' ',$helix[$i]);
	$start=$split[4].'_'.$split[3];
	$end=$split[7].'_'.$split[3];
	$helix{$start}=1;
	$helix{$end}=1;
	}
open(A,$ss_assign);
@data=<A>;
close(A);
for($i=0;$i<scalar@pp;$i++)
	{
	@split1=split(' ',$pp[$i]);
	$start=$split1[4]-1;
	$end=$split1[7]+1;
	$length=$end-$start+1;
	$chain=$split1[3];
	for($j=3;$j<scalar@data-13;$j++)
		{
		@split=split(' ',$data[$j]);
		if(((substr($data[$j],12,4)>=$start)&&(substr($data[$j],36,4)<=$end)&&(substr($data[$j],43,1) eq $chain)))
			{
			splice(@data,$j,1);$j--;
			}
		if(($split[7]>$end)&&($split[9] eq $chain))
			{
			last;
			}
		}
	}
push(@line,3);
for($j=0;$j<scalar@data;$j++)
	{
	if($data[$j]=~m/$protein/)
		{
		push(@line,$j);
		}
	}
for($k=0;$k<scalar@line-1;$k++)
	{
	$hi=0;@tab=();$hi1=0;@first=();
	$start=$line[$k];
	$end=$line[$k+1];
	for($m=$start+1;$m<$end;$m++)
		{
		@split=split(' ',$data[$m]);
		if($data[$m]=~m/P U U/)
			{
			substr($data[$m],0,1)='U';
			}
		if(((($split[13]>$t1)&&($split[13]<$t2))||(($split[-1]>$r1)&&($split[-1]<$r2)))&&($split[13]<240)&&(($split[14]>$h1)&&($split[14]<$h2)))
			{
			substr($data[$m],0,1)='S';
			$hi++;
			}
		push(@tab,$data[$m]);push(@first,substr($data[$m],0,1));
		}
	#print "$start\t$end\n";
	if($hi>1)
		{
		#print "$start\t$end\nori\n@tab\n";
		for($ch=(scalar@tab-1);$ch>=0;$ch--)
			{
			if($tab[$ch] =~m/U U U/)
				{
				splice(@tab,$ch,1);splice(@first,$ch,1);
				$ch=scalar@tab;
				}
			else
				{
				last;
				}
			}
		for($ch=0;$ch<scalar@tab;$ch++)
			{
			if($tab[$ch] =~m/U U U/)
				{
				splice(@tab,$ch,1);splice(@first,$ch,1);
				$ch=-1;
				}
			else
				{
				last;
				}
			}
		$l++;
		}
	else
		{
		$hi=0;@tab=();$hi1=0;@first=();
		}
	if(scalar@tab>0)
		{
		push(@tab,'U U U   61    69 G    70 A    71 R    72 P A   360.00  4.00  360.00   33.34  3.73');push(@first,'F');
		}
	for($ch=0;$ch<scalar@tab;$ch++)
		{
		@split=split(' ',$tab[$ch]);
		@split1=split(' ',$tab[$ch-1]);
		if(($tab[$ch]=~m/S U U/)&&($ch==0))
			{
			push(@tab1,$tab[$ch]);push(@twist,$split[13]);push(@first1,$first[$ch]);
			$number++;
			}
		elsif(($tab[$ch]=~m/S U U/)&&(($split[3]-$split1[3])==1)&&($ch>0))
			{
			push(@tab1,$tab[$ch]);push(@twist,$split[13]);push(@first1,$first[$ch]);
			$number++;
			}
		else
			{
			strand(@tab1);@tab1=();$number=0;@twist=();@first1=();
			}
		}
	}
@line=();
#print @strands;
#print "\n\n";
for($ch=0;$ch<scalar@strands;$ch++)
	{
	$rm1=0;$rm2=0;
	#print "$strands[$ch]\n";
	@split=split(' ',$strands[$ch]);
	$start=$split[4].'_'.$split[3];
	$end=$split[7].'_'.$split[3];

	$start1=($split[4]-2).'_'.$split[3];
	$end1=($split[7]-1).'_'.$split[3];
	@split1=split(' ',$nh{$start1});
	@split2=split(' ',$nh{$end1});
	#print "$start\t$helix{$start}\t$end\t$helix{$end}\n";
	#print "$nh{$start1}$nh{$end1}\n\n";
	if(($helix{$start}==1)||(($split1[10]>=77.9)&&($split1[10]<=114.9)&&($split[-1]>3)))
		{
		$rm1=1;
		if($split[-1]>3)
			{
			#print $strands[$ch];
			$split[4]=$split[4]+1;
			$split[-1]--;
			$reschain=$split[4].$split[3];
			$string=swrite(<<'END' ,'Strand',$split[1],$aminoacid{$reschain},$split[3],$split[4],$split[5],$split[6],$split[7],$split[8],$split[9]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$strands[$ch]=$string;
			}
		else
			{
	#		print "Removed1\t$ch\t$strands[$ch]";
			splice(@strands,$ch,1);splice(@remove,$ch,1);$ch--;
			next;
			}
		}
	if(($helix{$end}==1)||(($split2[10]>=77.9)&&($split2[10]<=114.9)&&($split[-1]>3)))
		{
		$rm2=1;
		if($split[-1]>3)
			{
			$split[7]=$split[7]-1;
			$split[-1]--;$reschain=$split[7].$split[3];
			$string=swrite(<<'END' ,'Strand',$split[1],$split[2],$split[3],$split[4],$aminoacid{$reschain},$split[6],$split[7],$split[8],$split[9]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$strands[$ch]=$string;

			}
		else
			{
	#		print "Removed2\t$ch\t$strands[$ch]";
			splice(@strands,$ch,1);splice(@remove,$ch,1);$ch--;
			next;
			}
		}
	}

for($ch=0;$ch<scalar@strands-1;$ch++)
	{
	@split=split(' ',$strands[$ch]);
	@split1=split(' ',$strands[$ch+1]);
	if(($split[7]==$split1[4])&&($split[3] eq $split1[3]))
		{
		$ch1=$ch+1;$turn[$ch].=$split[7].',';
		$split[7]=$split1[7];
		$split[-1]=$split[7]-$split[4]+1;
		$string=swrite(<<'END' ,'Strand',$split[1],$split[2],$split[3],$split[4],$split1[5],$split[6],$split[7],$split[8],$split[9]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
		$strands[$ch]=$string;splice(@strands,$ch1,1);splice(@remove,$ch1,1);$ch--;
		}
	}
chomp(@strands);
open(A,$ss_assign);
@data=<A>;
close(A);
open(C,">>$assignment");
for($ch=0;$ch<scalar@strands;$ch++)
	{
	@split=split(' ',$strands[$ch]);
	$ch1=$ch+1;
	$string=swrite(<<'END' ,'Strand',$ch1,$split[2],$split[3],$split[4],$split[5],$split[6],$split[7],$split[8],$split[9]);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
chomp($string);
	if($remove[$ch] eq 'yes')
		{
		splice(@strands,$ch,1);splice(@remove,$ch,1);$ch--;
		}
	else
		{
		if($turn[$ch] ne '')
			{
			print C "$string\tTwist in strand at $turn[$ch]\n";
			for($k=0;$k<scalar@data;$k++)
				{
				@split1=split(' ',$data[$k]);
				if(($split1[4]>=$split[4]-1)&&($split1[10]<=$split[7]+1)&&($split1[12] eq $split[3]))
					{
					substr($data[$k],0,5)='S U U';
					}
				if(($split1[10]>$split[7]+1)&&($split1[12] eq $split[3]))
					{
					last;
					}
				}
			}
		else
			{
			print C "$string\n";
			for($k=0;$k<scalar@data;$k++)
				{
				@split1=split(' ',$data[$k]);
				if(($split1[4]>=$split[4]-1)&&($split1[10]<=$split[7]+1)&&($split1[12] eq $split[3]))
					{
					substr($data[$k],0,5)='S U U';
					}
				if(($split1[10]>$split[7]+1)&&($split1[12] eq $split[3]))
					{
					last;
					}
				}
			}
		}
	}
close(C);
open(A,">$ss_assign");
print A @data;
close(A);
sub strand
	{
	local(@tab2)=@_;
	#print "prasun\n@tab2\n";
	if(scalar@tab2>1)
		{
		$twist1=sum(@twist)/@twist;
		$first1=join('',@first1);
		if($twist1> 180)
			{
			$k1++;
			@split1=split(' ',$tab2[0]);
			@split2=split(' ',$tab2[-1]);
			$start1=$split1[6];$ress=$split1[7];$rese=$split2[9];
			$end1=$split2[8];$chain=$split1[12];
			$res1=$ocode{$ress};$res2=$ocode{$rese};
			$length1=$end1-$start1+1;
			if($first1!~m/S/)
				{
				$remove[$s]='yes';
				}
			else
				{
				$remove[$s]='no';
				}
			$string = swrite(<<'END' ,'Strand',$k1,$res1,$chain,$start1,$res2,$chain,$end1,'0',$length1);
@<<<<<<<<<<@>>@>>>@>@>>>>@>>>>@>@>>>>@>>@>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
END
			$strands[$s]=$string;$s++;#print $string;
			}
		}
	$number=0;@tab2=();
	}
