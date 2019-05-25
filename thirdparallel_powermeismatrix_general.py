### aqui será mostrado como faremos o programa que trata de construir matrizes de uma forma geral, dada as geometrias de empilhamento que foram estudadas!!!
# -*- coding: utf-8 -*-

import numpy 
import threading
import math
import multiprocessing

from numpy import matrix,sqrt
from multiprocessing.pool import ThreadPool
from multiprocessing import Pool
from joblib import Parallel, delayed
from functools import partial

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################


# aqui será escolhida o ordenamento dos elementos na nano
## isto tem que ser revisto para fazer acordo com o programa PowerMEIS

def estruturas(alloytudo,cascaptpd,core,cs,dl):

	cnum=4.0

	if (alloytudo=='s' ):
		radiuspt=(float(input('Tamanho do Caroço de PtPd')))
		radiuspd=radiuspt
		radiusc=radiuspt+cs
		radiuscore=radiuspt
		radiusshell=radiuspt
		corenumber= 3
		shellnumber=3


#	print (radiusc,radiusshell,radiuscore)

	if(alloytudo == 'n'):


        	if(cascaptpd=='s' and core=='pt'):
			
        	    	radiuspt=float(input('Tamanho do Caroço de Pt'))
        	    	radiuspd=float(input('Tamanho da Casca de PtPd '))
        	    	radiusc=radiuspt+cs
        	    	radiuscore=radiuspt
        	    	radiusshell=radiuspd
        	    	corenumber=1
        	    	shellnumber=3.0


        	if(cascaptpd=='s' and core=='pd'):
        	    	radiuspd=float(input('Tamanho do Caroço de Pd '))
        	    	radiuspt=float(input('Tamanho da Casca de PtPd '))
        	    	radiusc=radiuspt+cs
        	    	radiuscore=radiuspd
        	    	radiusshell=radiuspt
        	    	corenumber=2
        	    	shellnumber=3.0




        	if(cascaptpd=='n'):

        	    	core=float(input('Tamanho do Caroço de pd ou pt?? '))

        	if(core=='pd'):
        	        radiuspd=float(input('Tamanho do Caroço de Pd '))
        	        radiuspt= float(input('Tamanho da Casca de Pt'))
        	        radiusc=radiuspt+cs
        	        radiuscore=radiuspd
        	        radiusshell=radiuspt
        	        corenumber=2
        	        shellnumber=1


        	if(core=='pt'):
        	        radiuspt=float(input('Tamanho do Caroço de Pt '))
        	        radiuspd=float(input('Tamanho da Casca de Pd'))
        	        radiusc=radiuspd+cs
        	        radiuscore=radiuspt
        	        radiusshell=radiuspd
        	        corenumber=1
        	        shellnumber=2

#	if(alloytudo == 'n' and cascaptpd == 'n' and (core =='pd' or core == 'pt')):
#		return [radiusc, radiuscore, radiusshell, corenumber, shellnumber, radiuspt, radiuspd,cnum]

#	if(alloytudo =='s'):
	return [radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum]



###########################################################################################################################

# fazer uma função que determine se a posição x,y,z está dentro do core-shell
# this function is to generate spherical nanoparticles

def incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z):

	div=float(radiusc/dl)


#	x[ii]= float((((float((1-(g-1)*sqrt(3.0)+Npiy)))*div-i+2)))
#	y[ii]= float(((float((-1+1-(g-1)+2*(w-1)+Npiy))*div-j)+2))
#	z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))
#	print (x[ii])

#	print(corenumber)


	aa=0

# core

	if ( ((x[ii]/(div))**2+(y[ii]/(div))**2+(z[ii]/(div))**2) <= (((float(radiuscore))/(float(radiusc)))**2) ):
					#			A[a][b]=corenumber
								index=corenumber
								aa=1
#								print(corenumber, "corenum")
# shell

	if ( (((x[ii]/(div))**2+(y[ii]/(div))**2+ (z[ii]/(div))**2) >= (((float(radiuscore))/(float(radiusc)))**2)) and (((x[ii]/(div))**2+(y[ii]/(div))**2+ (z[ii]/(div))**2) <= (((float(radiusshell))/(float(radiusc)))**2))):
#						                A[a][b]=shellnumber
								index=shellnumber
								aa=1

# carbon shell

	if ( (((x[ii]/((div)))**2+(y[ii]/((div)))**2+ (z[ii]/((div)))**2) >= (((float(radiusshell))/(float(radiusc)))**2)) and (((x[ii]/((div)))**2+(y[ii]/((div)))**2+ (z[ii]/((div)))**2) <= ((1.0)**2))):
#		        			        	A[a][b]=cnum
								index=cnum
								aa=1



	if(aa==0):
		index=0


	if(aa==1):
		pass
#		print (index,'index')


	return index,aa

###########################################################################################################################



# fazer uma função que determine se a posição x,y,z está dentro do core-shell
# this function is to generate interpenetrated spherical nanoparticles with a surfactant shell!!!!!!!!

def incoreshell_interpenesphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,rcarbonpt,radiuspt,rcarbonpd,radiuspd,distcentro,alloynumber):

	div=float(radiusc/dl)


	aa=0



# FIRST ELEMENT!

	if(radiuspt>radiuspd):
		largerr = rcarbonpt
	else:
		largerr = rcarbonpd


	xx1 = ((radiuspt + (largerr-radiuspt))/dl-i); 
#	xx1 = ((radiuspt + (rcarbonpt-radiuspt))/dl-i); 
	yy1 = ((radiuspt+(rcarbonpt-radiuspt))/dl-j); 
	zz1 = ((radiuspt+(largerr-radiuspt))/dl-k+1);

#	print(xx1)

# FIRST ELEMENT SHELL

	xx2 = ((radiuspt+(largerr-radiuspt))/dl-i) 
#	xx2 = ((radiuspt+(rcarbonpt-radiuspt))/dl-i) 
	yy2 = ((radiuspt+(rcarbonpt-radiuspt))/dl-j) 
	zz2 = ((radiuspt+(largerr-radiuspt))/dl-k+1)



# SECOND ELEMENT !

	xx3 = ((radiuspd+(largerr-radiuspd))/dl-i+1) 
#	xx3 = ((radiuspd+(rcarbonpd-radiuspd))/dl-i+1) 
	yy3 = ((radiuspd+(rcarbonpd-radiuspd)+distcentro)/dl-j+1)
	zz3 = ((radiuspd+(largerr-radiuspd))/dl-k+1)

# SECOND ELEMENT SHELL

	xx4 = ((radiuspd+(largerr-radiuspd))/dl-i+1)
#	xx4 = ((radiuspd+(rcarbonpd-radiuspd))/dl-i+1)
	yy4 = ((radiuspd+(rcarbonpd-radiuspd)+distcentro)/dl-j+1)
	zz4 = ((radiuspd+(largerr-radiuspd))/dl-k+1)

	em=1

	if (distcentro>rcarbonpt+rcarbonpd):


		if ( ((xx1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+ 
(zz1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2) <= (((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2) ):
			index=corenumber
			aa=1
#			print (aa)


		elif ((xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 >= ((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
and (xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 <= 1.0):
			index=cnum
			aa=1

		elif ((xx3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2):
			index=shellnumber
			aa=1

		elif ((xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 >= ((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 
and (xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 <= 1.0):
			index=cnum
			aa=1








	if (em == 1 and (xx1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2<=((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
 and (xx3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2>((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 
and (((xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2)+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2>1.0 or (xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 
+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<(((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2) ) ) :

		index=corenumber
		aa=1

		em=1
#xy(i,j+(k-1)*(Dy1))=1
 
	elif (em == 1 and (xx1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 
+(zz1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2<=((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2  
and ((xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2>=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 
and (xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 
+(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<=1.0)):
		em=1

		index=corenumber
		aa=1

#xy(i,j+(k-1)*(Dy1))=1



	elif (em == 1 and (xx1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 
+(zz1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2<=((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
and (xx3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 
+(zz3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2):
		em=1
		index=alloynumber
		aa=1


#xy(i,j+(k-1)*(Dy1))=3

# NOW FOR PT CARBON SHELL, i changed here

	if ((em == 1 and (xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2
+(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2>=((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
and (xx2/((radiuspt + (rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2<= 1.0 and (xx3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(yy3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2>=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 
and (((xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 )
or ((xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2)>=1.0))):
		em=1
		index= cnum
		aa=1


#xy(i,j+(k-1)*(Dy1))=4



	elif (em == 1 and (xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 >=((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
 and (xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2<= 1.0 
 and (xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 + 
(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2>=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 
and (xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<=1.0):
		em=1
		index=cnum
		aa=1
	

#xy(i,j+(k-1)*(Dy1))=4


	elif (em == 1 and (xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 >=((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
and (xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 <= 1.0 
and (xx3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<=(((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2) ):
		em=1
		index= shellnumber
		aa=1


#xy(i,j+(k-1)*(Dy1))=2



# just PD SPHERE! i changed here


	if (em == 1 and (xx3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz3/((radiuspd+(rcarbonpd-radiuspd))/dl))**2<=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 
 and (xx1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2>((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
and (((xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2< ((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
or ((xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2>1.0)))):

		em=1
		index=shellnumber
		aa=1

#xy(i,j+(k-1)*(Dy1))=2



# JUST PD CARBON SHELL, i changed here

	if ((em == 1 and (xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 >=((radiuspd)/(radiuspd+(rcarbonpd-radiuspd)))**2 
and (xx4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+(yy4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2+
(zz4/((radiuspd+(rcarbonpd-radiuspd))/dl))**2 <= 1.0 and (xx1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(yy1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz1/((radiuspt+(rcarbonpt-radiuspt))/dl))**2>((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 
and ( ((xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
(zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2 < ((radiuspt)/(radiuspt+(rcarbonpt-radiuspt)))**2 ) 
or (( (xx2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+(yy2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2+
((zz2/((radiuspt+(rcarbonpt-radiuspt))/dl))**2) > 1.0))))):
	
		em=1	
		index=cnum
		aa=1


#xy(i,j+(k-1)*(Dy1))=4






	if(aa==0):
		index=0


	if(aa==1):
		pass
#		print (index,'index')


	return index,aa




	return index,aa

###########################################################################################################################




###########################################################################################################################

# fazer uma função que determine se a posição x,y,z está dentro do core-shell
# this function is to generate spherical nanoparticles

def incoreshell_ellipsoid(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):

	div1=float(radiusxcs/dl)
	div2=float(radiusycs/dl)
	div3=float(radiuszcs/dl)


#	x[ii]= float((((float((1-(g-1)*sqrt(3.0)+Npiy)))*div-i+2)))
#	y[ii]= float(((float((-1+1-(g-1)+2*(w-1)+Npiy))*div-j)+2))
#	z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))
#	print (x[ii])

	aa=0




	if ( ((x[ii]/(div1))**2+(y[ii]/(div2))**2+(z[ii]/(div3))**2) <= (((float(radiuscore))/(float(radiusc)))**2) ):
					#			A[a][b]=corenumber
								index=corenumber
								aa=1
# shell

	if ( (((x[ii]/(div1))**2+(y[ii]/(div2))**2+ (z[ii]/(div3))**2) >= (((float(radiuscore))/(float(radiusc)))**2)) and (((x[ii]/(div1))**2+(y[ii]/(div2))**2+ (z[ii]/(div3))**2) <= (((float(radiusshell))/(float(radiusc)))**2))):
#						                A[a][b]=shellnumber
								index=shellnumber
								aa=1

# carbon shell

	if ( (((x[ii]/((div1)))**2+(y[ii]/((div2)))**2+ (z[ii]/((div3)))**2) >= (((float(radiusshell))/(float(radiusc)))**2)) and (((x[ii]/((div1)))**2+(y[ii]/((div2)))**2+ (z[ii]/((div3)))**2) <= ((1.0)**2))):
#		        			        	A[a][b]=cnum
								index=cnum
								aa=1






	if(aa==0):
		index=0


	if(aa==1):
		pass
#		print (index,'index')


	return index,aa

###########################################################################################################################




###########################################################################################################################

# fazer uma função que determine se a posição x,y,z está dentro do core-shell
# this function is to generate cubic nanoparticles

def incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z):

	div=float(radiusc/dl)


#	print (x[ii])

	aa=0

# core


	acore=(x[ii]/(div))**2 <= (((float(radiuscore))/(float(radiusc)))**2) and (y[ii]/(div))**2 <= (((float(radiuscore))/(float(radiusc)))**2) and (z[ii]/(div))**2<= (((float(radiuscore))/(float(radiusc)))**2)
	
	if ( acore ):
					#			A[a][b]=corenumber
								index=corenumber
								aa=1
# shell


#	a1x= (x[ii]/(div))**2>= ((float(radiuscore))/(float(radiusc)))**2 and (x[ii]/(div))**2<= ((float(radiusshell))/(float(radiusc)))**2



#	a1y= (y[ii]/(div))**2 >= ((float(radiuscore))/(float(radiusc)))**2 and  (y[ii]/(div))**2 <= ((float(radiusshell))/(float(radiusc)))**2 


#	a1z= (z[ii]/(div))**2 >= ((float(radiuscore))/(float(radiusc)))**2  and  (z[ii]/(div))**2 <= ((float(radiusshell))/(float(radiusc)))**2


	ashell= (x[ii]/(div))**2<= ((float(radiusshell))/(float(radiusc)))**2 and  (y[ii]/(div))**2 <= ((float(radiusshell))/(float(radiusc)))**2  and  (z[ii]/(div))**2 <= ((float(radiusshell))/(float(radiusc)))**2

	if ( ashell and acore==False):
#	if ( 1==1):
#						                A[a][b]=shellnumber
								index=shellnumber
								aa=1

# carbon shell


#	a1= (x[ii]/(div))**2>= ((float(radiusshell))/(float(radiusc)))**2 and (y[ii]/(div))**2 >= ((float(radiusshell))/(float(radiusc)))**2 and (z[ii]/(div))**2 >= ((float(radiusshell))/(float(radiusc)))**2


#	a2= (x[ii]/(div))**2<=1.0 and (y[ii]/(div))**2 <= 1.0 and (z[ii]/(div))**2 <= 1.0	


	acarbon=(x[ii]/(div))**2<=1.0 and (y[ii]/(div))**2 <= 1.0 and (z[ii]/(div))**2 <= 1.0

	if ( ashell==False and acarbon==True):
#		        			        	A[a][b]=cnum
								index=cnum
								aa=1



	if(aa==0):
		index=0


	if(aa==1):
		pass
#		print (index,'index')



	return index,aa

###########################################################################################################################



###########################################################################################################################

# fazer uma função que determine se a posição x,y,z está dentro do core-shell
# this function is to generate octahedrical nanoparticles

def incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z):

	div=float(radiusc/dl)


#	print (x[ii])

	aa=0

# core

	acore=abs(x[ii]/(div))+abs(y[ii]/(div))+abs(z[ii]/(div))<= (((float(radiuscore))/(float(radiusc)))**2)


	if ( acore ):

					#			A[a][b]=corenumber
								index=corenumber
								aa=1
# shell


	ashell= abs(x[ii]/(div))+abs(y[ii]/(div))+abs(z[ii]/(div))<= (((float(radiusshell))/(float(radiusc)))**2)

	if ( ashell and acore==False ):

#						                A[a][b]=shellnumber
								index=shellnumber
								aa=1

# carbon shell

	acarbon= abs(x[ii]/(div))+abs(y[ii]/(div))+abs(z[ii]/(div))<= (((float(radiusc))/(float(radiusc)))**2) 


	if ( acarbon and ashell == False ):

#		        			        	A[a][b]=cnum
								index=cnum
								aa=1


	if(aa==0):
		index=0


	if(aa==1):
		pass
#		print (index,'index')


	return index,aa

###########################################################################################################################

###########################################################################################################################

# fazer uma função que determine se a posição x,y,z está dentro do core-shell
# this function is to generate octahedrical nanoparticles

def incoreshell_triplate(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,base,width,height,base2,width2,cs):


#### here are described the equations to form a triangular plate nanoparticle
### it must be repaired

	div=float(radiusc/dl)


	basenew= base2+cs
	widthnew=width2+cs

	aa=0

	acore1=False
	acore2=False


	if( y[ii]/div <= width/(2*widthnew)+abs(x[ii]/div)*(width*(basenew)/base*(widthnew)) and y[ii]>=0):
#		print (x[ii],base/2)

		acore1=True


#	print (x[ii],base/2)


#    	        		x[ii]= float((((float(2.0*(g-1)+1+(h-1)))*float(radiusc)/float(dl)-float(i)

#	if( (x[ii]+float(i))/div <=x[ii]/div<= base/(2*basenew)):

#		if( (y[ii]+float(j))/div <=width/widthnew-i*(2*width*(basenew)/(base*widthnew))):
#			print (x[ii],base/2)

#			acore1=True

#	if( -base/(2*basenew) <= x[ii]/div <= 0 ):

#		if( y[ii]/div <=width/widthnew+(x[ii]/div)*(2*width*(basenew)/(base*widthnew))):

#			acore2= True
#


	acore2=False

	if ( acore1 or acore2 ):

#								print (x[ii],base/2)


					#			A[a][b]=corenumber
								index=corenumber
								index=0
								aa=1

#								print(index)



	ashell1=False
	ashell2=False



	if( 0 <=x[ii]<= base2/2):

		if( y[ii]<=width2-x[ii]*(2*width2/base2)):

			ashell1=True

	if( -base2/2 <= x[ii] <= 0 ):

		if( y[ii]<=width2+x[ii]*(2*width2/base2)):

			ashell2= True




	aant=True
	
	if (acore1==False and acore2==False):

		aant=False

	if ( ashell1 or ashell2 and ( aant==False) ):


					#			A[a][b]=corenumber
								index=shellnumber
								aa=1
#								print(index)

	acnum1=False
	acnum2=False

	basenew= base2+cs
	widthnew=width2+cs


	if( 0 <=x[ii]<= basenew/2):

		if( y[ii]<=widthnew-x[ii]*(2*widthnew/basenew)):

			acnum1=True

	if( -base2/2 <= x[ii] <= 0 ):


		if( y[ii]<=widthnew+x[ii]*(2*widthnew/basenew)):
			acnum2= True

	aant=True
	
	if (ashell1==False and ashell2==False):

		aant=False

	if ( acnum1 or acnum2 and ( aant==False) ):



					#			A[a][b]=corenumber
								index=cnum
								index=0
								aa=1

#								print(index)



# shell



#	print (x[ii])


# core


	if(aa==0):
		index=0


	if(aa==1):
		pass
#		print (index,'index')

#	print(index)
	return index,aa

###########################################################################################################################



# condições iniciais dos arrays criados

def zeroarray(Npix,Npiy,Ntam):

#jj=1

#w=Npiy
#g=Npix


	if(Npix>Npiy):
		Nmaior=Npix
	else:
		Nmaior=Npiy



	Npart=0
	for h in range(0,(Nmaior+1)):
		Npart=Npart+Nmaior**2


	Npartz=0
	for h in range(0,(Ntam+1)):
		Npartz=Ntam+h**2
	

	if(Npart> Npartz):
		pass

	else:
		Npart=Npartz

#	print Npart,'Npart'

	x= [0 for i in range(Npart+1)]
	y= [0 for i in range(Npart+1)]
	z= [0 for i in range(Npart+1)]

	return [x, y, z, Npart] #=[x,y,z,Npart]

#[x, y, z, Npart] = zeroarray(10, 1)

###########################################################################################################################



### cria matriz e zera ela

def zeromatrix(Dx1,Dy1,Dz1):

	mult=(Dy1)*(Dz1)
#print (Dx1)

	A= [[0 for i in range(mult+1)] for j in range(Dx1+1)]


#print A[Dx1][mult]
#print len(A)



#A=[[]]


	for k in range(1,(Dz1+1)):
    		for j in range(1,(Dy1+1)):
        		for i in range(1,Dx1+1):

 		        	A[i][j+(k-1)*Dy1]=0
#            A[i][j]=0
#        print i


	return A
###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções



def jonderstruct1 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs):




	w=0
	g=0
	h=0
#	print(corenumber,'aqui')

	for k in range(zini,zfin+1):
		for j in range(yini,yfin+1):
		        for i in range(xini,xfin+1):
		        	a=int(i)
		        	b=int(j+(k-1)*Dy1)
		        	ii=1

#############

		        	xx1=((1)*float(radiusc/dl)-i+1)
		        	yy1= ((1)*float(radiusc/dl)-j+1)
		        	zz1=(float(radiusc/dl)-k+1)
    
		        	xx2=((1.0)*float(radiusc/dl)-i+1)
		        	yy2= ((3.0)*float(radiusc/dl)-j+1)
		        	zz2=(float(radiusc/dl)-k+1)

		        	xx3=((3.0)*float(radiusc/dl)-i+1)
		        	yy3= ((1.0)*float(radiusc/dl)-j+1)
		        	zz3=(float(radiusc/dl)-k+1)

		        	xx4=((3.0)*float(radiusc/dl)-i+1)
		        	yy4= ((3.0)*float(radiusc/dl)-j+1)
		        	zz4=(float(radiusc/dl)-k+1)

#!!! second layer, one sphere

		        	xx5=((2.0)*float(radiusc/dl)-i+1)
		        	yy5= ((2.0)*float(radiusc/dl)-j+1)
		        	zz5=((sqrt(2.0)+1.0)*float(radiusc/dl)-k+1)

#!!! third layer, four spheres


		        	xx6=(1*float(radiusc/dl)-i+1)
		        	yy6= (1*float(radiusc/dl)-j+1)
		        	zz6=((sqrt(3.0)+2.0)*float(radiusc/dl)-k+1)

		        	xx7=(1*float(radiusc/dl)-i+1)
		        	yy7= (3*float(radiusc/dl)-j+1)
		        	zz7=((sqrt(3.0)+2.0)*float(radiusc/dl)-k+1)

		        	xx8=(3*float(radiusc/dl)-i+1)
		        	yy8= (1*float(radiusc/dl)-j+1)
		        	zz8=((sqrt(3.0)+2.0)*float(radiusc/dl)-k+1)
    
		        	xx9=(3*float(radiusc/dl)-i+1)
		        	yy9= (3*float(radiusc/dl)-j+1)
		        	zz9=((sqrt(3.0)+2.0)*float(radiusc/dl)-k+1)


		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )


#            print radiuscore,radiusshell,radiusc
            
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
			        	A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum



		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum
	

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum


		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum


		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum






		        	ii=ii+1






	return A


###########################################################################################################################


###########################################################################################################################

def jonderstruct2 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs):




	w=0
	g=0
	h=0
#	print(corenumber,'aqui')

	for k in range(zini,zfin+1):
		for j in range(yini,yfin+1):
		        for i in range(xini,xfin+1):
		        	a=int(i)
		        	b=int(j+(k-1)*Dy1)
		        	ii=1

		        	xx1=((1.0)*float(radiusc/dl)-i+1)
		        	yy1= ((1.0)*float(radiusc/dl)-j+1)
		        	zz1=(float(radiusc/dl)-k+1)
    

		        	xx2=((2.0+sqrt(3.0))*float(radiusc/dl)-i+1)
		        	yy2= ((1.0)*float(radiusc/dl)-j+1)
		        	zz2=(float(radiusc/dl)-k+1)



		        	xx3=((1.0)*float(radiusc/dl)-i+1)
		        	yy3= ((2.0+sqrt(3.0)) *float(radiusc/dl)-j+1)
		        	zz3=(float(radiusc/dl)-k+1)




		        	xx4=((2+sqrt(3.0))*float(radiusc/dl)-i+1)
		        	yy4= ((2+sqrt(3.0))*float(radiusc/dl)-j+1)
		        	zz4=(float(radiusc/dl)-k+1)




		        	xx5=((1.0+sqrt(2.0))*float(radiusc/dl)-i+1)
		        	yy5= ((1.0+sqrt(2.0))*float(radiusc/dl)-j+1)
		        	zz5=(float(radiusc/dl)-k+1)

#!!!!

		        	xx6=(1.0*float(radiusc/dl)-i+1)
		        	yy6= (1.0*float(radiusc/dl)-j+1)
		        	zz6=((2+sqrt(3.0))*float(radiusc/dl)-k+1)




		        	xx7=((2+sqrt(3.0))*float(radiusc/dl)-i+1)
		        	yy7= (1*float(radiusc/dl)-j+1)
		        	zz7=((2+sqrt(3.0))*float(radiusc/dl)-k+1)






		        	xx8=((1+sqrt(2.0))*float(radiusc/dl)-i+1)
		        	yy8= ((1)*float(radiusc/dl)-j+1)
		        	zz8=((1+sqrt(2.0))*float(radiusc/dl)-k+1)
    
		        	xx9=((1.0)*float(radiusc/dl)-i+1)
		        	yy9= ((2+sqrt(3.0))*float(radiusc/dl)-j+1)
		        	zz9=((2+sqrt(3.0))*float(radiusc/dl)-k+1)



		        	xx10=((2+sqrt(3.0))*float(radiusc/dl)-i+1)
		        	yy10= ((2+sqrt(3.0))*float(radiusc/dl)-j+1)
		        	zz10=((2+sqrt(3.0))*float(radiusc/dl)-k+1)



		        	xx11=((1.0+sqrt(2.0))*float(radiusc/dl)-i+1)
		        	yy11= ((2.0+sqrt(3.0))*float(radiusc/dl)-j+1)
		        	zz11=((1.0+sqrt(2.0))*float(radiusc/dl)-k+1)

		        	xx12=((2.0+sqrt(3.0))*float(radiusc/dl)-i+1)
		        	yy12= ((1.0+sqrt(2.0))*float(radiusc/dl)-j+1)
		        	zz12=((1.0+sqrt(2.0))*float(radiusc/dl)-k+1)

		        	xx13=((1.0+sqrt(2.0))*float(radiusc/dl)-i+1)
		        	yy13= ((1.0+sqrt(2.0))*float(radiusc/dl)-j+1)
		        	zz13=((2.0+sqrt(3.0))*float(radiusc/dl)-k+1)





		        	xx14=(1*float(radiusc/dl)-i+1)
		        	yy14= ((1.0+sqrt(2.0))*float(radiusc/dl)-j+1)
		        	zz14=((1.0+sqrt(2.0))*float(radiusc/dl)-k+1)






		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )


#            print radiuscore,radiusshell,radiusc
            
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
			        	A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum



		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum
	

		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum


		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum


		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum



		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum




		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum




		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum





		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum




		        	div=float(radiusc)/float(dl)
		        	cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )

	
		        	if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
		        		A[a][b]=corenumber

		        	if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
		        		A[a][b]=shellnumber

		        	if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
		        		A[a][b]=cnum








		        	ii=ii+1






	return A


###########################################################################################################################

def para_rec(j,k,xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs):


#        print(radiuscore, radiusc, "dimensions")


#        print(j,"j")

        iii=0
        jjj=0
        ind=0
        listind=[]

        for i in range(xini,xfin+1):
        	a=int(i)
        	b=int(j+(k-1)*Dy1)
        	ii=1
	        for h in range(zini,Ntam+1):

		        for g in range(zini,Npiy+1):

		
#					        w_queue= Queue()
#					        w_queue= Queue()
			    	
			        for w in range(zini,Npix+1):
#				print a,b

#    xx1=(basalradius/dl-i+1)
#   yy1= (basalradius/dl-j+1)
#    zz1=(perpradius/dl-k+1)


	    	        		x[ii]= float(((float(radiusc*(2*(w-1)+1))/float(dl)-float(i)+2.0)))
	    	        		y[ii]= float(((float(radiusc*(2*(g-1)+1))/float(dl)-float(j))+2.0))
	    	        		z[ii]= float(((float(radiusc*(2*(h-1)+1))/float(dl)-float(k)+2.0)))
#	    	        		print(x[ii],ii,len(x))


#	    	        		print(i,j,k, "i,j,k")



#	    	        		print(w,g,h,"w,g,h")






	    	        		if (shape=='sphere'):

	    	        			[index,aa] = incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)
	
	    	        			if (aa==1):
#	    	        				A[i][j+(k-1)*Dy1] = index
	    	        				iii=i
	    	        				jjj=j+(k-1)*Dy1
#	    	        				if(index==2):
#	    	        					print("indice igual a 1")
	    	        				ind=index
	    	        				listind.append([iii,jjj,ind])
#	    	        				print(ind, "eh o indice")



#	    	        				print(index)



	    	        		if (shape=='cubic'):
	
	    	        			[index,aa] = incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

	    	        			if (aa==1):
	    	        				A[i][j+(k-1)*Dy1] = index



	    	        		if (shape=='octahedron'):	

	    	        			[index,aa] = incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

	    	        			if (aa==1):
		    	        			A[i][j+(k-1)*Dy1] = index



	    	        		if (shape=='triplate'):	


	    	        			[index,aa] = incoreshell_triplate(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,base,width,height,base2,width2,cs)

	    	        			if (aa==1):
	    	        				A[i][j+(k-1)*Dy1] = index


	    	        		ii=ii+1


#        return A
#        if(ind==1):
#        	print(ind)
#        print(listind)

#        return iii,jjj,ind
	
        return listind


###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções

def rectangular (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs):

############################################################################################################################

## TRYING TO PARALLELIZE THIS PART!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

############################################################################################################################

#	nThreads = 4


#	pool = ThreadPool(processes=nThreads)


	num_cores = multiprocessing.cpu_count()
#	print(num_cores, "cores")

#	w=0
#	g=0
#	print(corenumber,'aqui')
#	print(1)

	for k in range(zini,zfin+1):

#		print(zini, zfin+1, "zs")


#		for j in range(yini,yfin+1):
#		inputs = range(yini,yfin+1)
#		parts= partial(para_rec,k=k,xini=xini,xfin=xfin,yini=yini,yfin=yfin,zini=zini,zfin=zfin,Npix=Npix,Npiy=Npiy,Ntam=Ntam,Dx1=Dx1,Dy1=Dy1,Dz1=Dz1,radiusc=radiusc,radiuscore=radiuscore,radiusshell=radiusshell,A=A,corenumber=corenumber,shellnumber=shellnumber,cnum=cnum,dl=dl,x=x,y=y,z=z,shape=shape,base=base,width=width,height=height,base2=base2,width2=width2,radiusx1=radiusx1,radiusy1=radiusy1,radiusz1=radiusz1,radiusx2=radiusx2,radiusy2=radiusy2,radiusz2=radiusz2,radiusxcs=radiusxcs,radiusycs=radiusycs,radiuszcs=radiuszcs,cs=cs)

		nThreads = 2*num_cores
		pool = ThreadPool(processes=nThreads)

		threadsrange=int((yfin+1)/nThreads)
#		print(threadsrange)



		proc = []



		for kkk in range(nThreads):					        	
			for j in range(kkk*threadsrange,(kkk+1)*threadsrange):
#				print(j,"eh o j")	

				proc.append(pool.apply_async(para_rec, (j,k,xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs)))

			
#		print(len(proc), "length of proc")

		for j in range(len(proc)):
			proc[j].wait()
			


		out = []

		for j in range(len(proc)):
			out.append(proc[j].get())


#		print(len(out))
#		print()


#		print(len(out), "comp1",len(out[0]), "comp2")
		

#		w=out[0]
#		print(w[0])
#		print(len(out),"length of out")
#		print(len(w),"length of w")



		for j in range(len(out)):
			ww=out[j]
#			print(w)
			if(ww!=[]): 
#				print(w)
				for o in range(len(ww)):
					vv=ww[o]
					for dd in range(len(vv)):
						A[vv[0]][vv[1]]=vv[2]				


#				xxxx=1
#			print(len(w))
#			print(len(out))
#			print (w[2])
#			for o in range(len(w)):
#			if(w[2]!=0):
#			A[w[0]][w[1]]=w[2]				
#					if (w[2]==1):
#						print("w[2] igual a 1",w[0],w[1])




#		A = Parallel(n_jobs=num_cores)(delayed(parts)(j) for j in inputs)


		
#		para_rec (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs)



#	print(A[21][14],"esse")

	return A


###########################################################################################################################




###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções

def columnpiling (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs):




	w=0
	g=0
#	print(corenumber,'aqui')

	for k in range(zini,zfin+1):
		for j in range(yini,yfin+1):
		        for i in range(xini,xfin+1):
		        	a=int(i)
		        	b=int(j+(k-1)*Dy1)
		        	ii=1
			        for h in range(zini,Ntam+1):

			    	
#				print a,b

#    xx1=(basalradius/dl-i+1)
#   yy1= (basalradius/dl-j+1)
#    zz1=(perpradius/dl-k+1)

	    	        		x[ii]= float(((float(radiusc)/float(dl)-float(i)+2.0)))
	    	        		y[ii]= float(((float(radiusc)/float(dl)-float(j))+2.0))
	    	        		z[ii]= float(((float(radiusc*(2*(h-1)+1))/float(dl)-float(k)+2.0)))
#	    	        		print(x[ii],ii,len(x))

	    	        		if (shape=='sphere'):

	    	        			[index,aa] = incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

	    	        			if (aa==1):
		    	        			A[i][j+(k-1)*Dy1] = index



	    	        		if (shape=='cubic'):

	    	        			[index,aa] = incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

	    	        			if (aa==1):
		    	        			A[i][j+(k-1)*Dy1] = index



	    	        		if (shape=='octahedron'):	

	    	        			[index,aa] = incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

	    	        			if (aa==1):
		    	        			A[i][j+(k-1)*Dy1] = index



	    	        		if (shape=='triplate'):	


	    	        			[index,aa] = incoreshell_triplate(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,base,width,height,base2,width2,cs)

	    	        			if (aa==1):
	    	        				A[i][j+(k-1)*Dy1] = index


	    	        		ii=ii+1




	return A


###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções


def parallelpygen(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,xini,xfin,yini,yfin,zini,zfin,Dx1,Dy1,Dz1):

	mult=(Dy1)*(d_k)


	aux= [[0 for i in range(mult+1)] for j in range(Dx1+1)]

#print (Dx1)

#	A= [[0 for i in range(mult+1)] for j in range(Dx1+1)]








	for j in range(yini,yfin+1):

		for i in range(xini,xfin+1):

			a=int(i)

			b=int(j+(k-1)*Dy1)
			    	
#				print a,b
			ii=1

			for h in range(zini,Ntam+1):

				Npiyw=Npiy-(h-1)
				w=Npiyw
#					print h					
					        


				while (w > 0):

					Npixg=Npix-(h-1)
						
					g=Npixg

					while(g>0):

						x[ii]= float((((float(2.0*(g-1)+1+(h-1)))*float(radiusc)/float(dl)-float(i)+2.0)))
						y[ii]= float((((float(2.0*(w-1)+1+(h-1)))*float(radiusc)/float(dl)-float(j))+2.0))
						z[ii]= float((((float(sqrt(3.0)*(h-1)+1))*float(radiusc)/float(dl)-float(k)+2.0)))
	

						if (shape=='sphere'):

							[index,aa] = incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

							if (aa==1):
								aux[i][j+(k-1)*Dy1] = index



						if (shape=='cubic'):

							[index,aa] = incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

							if (aa==1):
								aux[i][j+(k-1)*Dy1] = index



						if (shape=='octahedron'):	

							[index,aa] = incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

							if (aa==1):
								aux[i][j+(k-1)*Dy1] = index



						if (shape=='triplate'):	

							[index,aa] = incoreshell_triplate(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,base,width,height,base2,width2,cs)

							if (aa==1):
								aux[i][j+(k-1)*Dy1] = index




						if (shape=='interpene'):	

							[index,aa] = incoreshell_interpenesphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,rcarbonpt,radiuspt,rcarbonpd,radiuspd,distcentro,alloynumber)

							if (aa==1):
								aux[i][j+(k-1)*Dy1] = index


#				    	        	while (g>0):
#								incoreshell(Npiy,Npix,w,g,i,j,k,ii,radiusc,radiuscore,radiusshell,A,a,b,corenumber,shellnumber,cnum)
#								!!!

#				    	        		div=float(radiusc/dl)


#				    	        		x[ii]= float((((float((1-(g-1)*sqrt(3.0)+Npiy)))*div-i+2)))
#				    	        		y[ii]= float(((float((-1+1-(g-1)+2*(w-1)+Npiy))*div-j)+2))
#				    	        		z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))


#								!!!
#				    	        		if(index != 0):

#					    	        		print(index)

						g=g-1
			            			
#							print ii,A[a][b]
						ii=ii+1
#							print g
					w=w-1
#





	return aux




def pygeneral (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs,rcarbonpt,radiuspt,rcarbonpd,radiuspd,distcentro,alloynumber):

	proc=[]
	n_threads=2
	i_thread=range(n_threads)

	pool = ThreadPool(processes=n_threads)


	print(corenumber,'aqui')
#	print(1)

	d_k=(zfin-zini)/n_threads

	for i_th in i_thread:
		for k in range(zini+i_th*d_k,zini+(+1+i_th)*d_k):
	
			proc.append(pool.apply_async(parallelpygen,(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,xini,xfin,yini,yfin,zini,zfin,Dx1,Dy1,Dz1)))

	for i in range(nThreads):
		proc.wait()

	matrix_out = []

	for i in range(nThreads):
		matrix_out.append(proc[i].get())

	for i_th in i_thread:
		for ll in range(zini+i_th*d_k,zini+(+1+i_th)*d_k):
			for j in range(yini,yfin+1):
				for i in range(xini,xfin+1):

					# there must be a formula for k!!!!
					k= (zini+(+1+i_th)*d_k)*i_th

					A[i][j+(k-1)*Dy1]=proc[i][j+(ll-1)*Dy1]



#print A[a][b],a,b

	return A

###########################################################################################################################


###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções
# descreve uma piramide com base hexagonal arbitrária


def hexpyramid (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


#	print(corenumber,'aqui')


	for v in range(1,5):
		for k in range(zini,zfin+1):
			for j in range(yini,yfin+1):
			        for i in range(xini,xfin+1):
			        	a=int(i)
			        	b=int(j+(k-1)*Dy1)
			        	ii=1
			        	for h in range (1, Ntam+1):	
					
			        		Npiyw=Npiy
	#	w=Npiyw/2

			        		Npixg=Npix
#	g=Npixg/2
			        		box=0


			        		if (h>1):


#			        		if( h> 1):

			        			Npiyw=(Npiy-(h-1))*2-1

			        		if(Npiyw<=0):
			        			Npiyw=Npiy+2

	
			        		for w in range (1,int(Npiyw/2+1)):

			        			box=box+1
			        			Npixg=Npixg-1

			        			if( h> 1):
			        				Npixg=2*(Npix-(h-1))
		
			        			if(Npixg<=0):
			        				Npixg=Npix+2


	
			        			for g in range(1,int(Npixg/2+1)):



			        				if (v==1):
			
			        					div=float(radiusc/dl)

			        					x[ii]= float((((float((1+(g-1)*sqrt(3.0)+Npiyw)))*div-i+2)))
			        					y[ii]= float(((float((+1+(g-1)+2*(w-1)))*div-j)+2))
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))


			        				if (v==2):

			        					div=float(radiusc/dl)

			        					x[ii]= float((((float((1-(g-1)*sqrt(3.0)+Npiyw)))*div-i+2)))
			        					y[ii]= float(((float((+1+(g-1)+2*(w-1)))*div-j)+2))
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))


	
			        				if (v==3 ):
		
			        					div=float(radiusc/dl)
	
		
			        					x[ii]= float((((float((1+(g-1)*sqrt(3.0)+Npiyw)))*div-i+2)))
			        					y[ii]= float(((float((-1+1-(g-1)+2*(w-1)+Npiyw))*div-j)+2))	
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))


			        				if (v==4 ):

	
			        					div=float(radiusc/dl)


			        					x[ii]= float((((float((1-(g-1)*sqrt(3.0)+Npiyw)))*div-i+2)))
			        					y[ii]= float(((float((-1+1-(g-1)+2*(w-1)+Npiyw))*div-j)+2))
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))



			        				if (shape=='sphere'):

			        					[index,aa] = incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

			        					if (aa==1):
			        						A[i][j+(k-1)*Dy1] = index



			        				if (shape=='cubic'):

			        					[index,aa] = incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

			        					if (aa==1):
			        						A[i][j+(k-1)*Dy1] = index



			        				if (shape=='octahedron'):	

			        					index = incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

			        					if (aa==1):
			        						A[i][j+(k-1)*Dy1] = index




	
			        				g=g-1




			        				ii=ii+1
			        			w=w-1


#print A[a][b],a,b

	return A

###########################################################################################################################


###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções

def hexmono (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


#	print(corenumber,'aqui')

	for v in range(1,5):
		for k in range(zini,zfin+1):
			for j in range(yini,yfin+1):
			        for i in range(xini,xfin+1):	
			        	a=int(i)
			        	b=int(j+(k-1)*Dy1)
			        	ii=1
			        	for h in range (1, Ntam+1):

			        		Npiyw=Npiy
#	w=Npiyw/2

			        		Npixg=Npix
#	g=Npixg/2
			        		box=0
						
#			        		print (Npiy)



#			        		print (int(Npix/2)+1)
			        		for w in range (1,int(Npiy/2)+1):

			        			box=box+1
			        			Npixg=Npixg-1
	

			        			for g in range(1,int(Npix/2+1)):



			        				if (v==1):
			
			        					div=float(radiusc/dl)


#	x(ii)= ((((1+(g-1)*Sqrt(3.0)+Npiy))*radiusc/dl-i+2))
#	y(ii)= (((+1+(g-1)+2*(w-1))*radiusc/dl-j)+2)
#	z(ii)= (((Sqrt(3.0)*(h-1)+1)*radiusc/dl-k+2))





			        					x[ii]= float((((float((1+(g-1)*sqrt(3.0)+Npiy)))*div-i+2)))
			        					y[ii]= float(((float((+1+(g-1)+2*(w-1)))*div-j)+2))
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))


			        				if (v==2):

			        					div=float(radiusc/dl)



#	x(ii)= ((((1-(g-1)*Sqrt(3.0)+Npiy))*radiusc/dl-i+2))
#	y(ii)= (((+1+(g-1)+2*(w-1))*radiusc/dl-j)+2)
#	z(ii)= (((Sqrt(3.0)*(h-1)+1)*radiusc/dl-k+2))

			        					x[ii]= float((((float((1-(g-1)*sqrt(3.0)+Npiy)))*div-i+2)))
			        					y[ii]= float(((float((+1+(g-1)+2*(w-1)))*div-j)+2))
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))

	
		
			        				if (v==3 ):
	
			        					div=float(radiusc/dl)

	


#	x(ii)= ((((1+(g-1)*Sqrt(3.0)+Npiy))*radiusc/dl-i+2))
#	y(ii)= (((-1+1-(g-1)+2*(w-1)+Npiy)*radiusc/dl-j)+2)
#	z(ii)= (((Sqrt(3.0)*(h-1)+1)*radiusc/dl-k+2))

			        					x[ii]= float((((float((1+(g-1)*sqrt(3.0)+Npiy)))*div-i+2)))
			        					y[ii]= float(((float((-1+1-(g-1)+2*(w-1)+Npiy))*div-j)+2))	
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))
	

			        				if (v==4 ):

	
			        					div=float(radiusc/dl)

#	x(ii)= ((((1-(g-1)*Sqrt(3.0)+Npiy))*radiusc/dl-i+2))
#	y(ii)= (((-1+1-(g-1)+2*(w-1)+Npiy)*radiusc/dl-j)+2)
#	z(ii)= (((Sqrt(3.0)*(h-1)+1)*radiusc/dl-k+2))



			        					x[ii]= float((((float((1-(g-1)*sqrt(3.0)+Npiy)))*div-i+2)))
			        					y[ii]= float(((float((-1+1-(g-1)+2*(w-1)+Npiy))*div-j)+2))
			        					z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))



			        				if (shape=='sphere'):

				    	        			[index,aa] = incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

				    	        			if (aa==1):
					    	        			A[i][j+(k-1)*Dy1] = index



			        				if (shape=='cubic'):

				    	        			[index,aa] = incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

				    	        			if (aa==1):
					    	        			A[i][j+(k-1)*Dy1] = index



			        				if (shape=='octahedron'):	

				    	        			index = incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

				    	        			if (aa==1):
					    	        			A[i][j+(k-1)*Dy1] = index



	
#			        				g=g-1



	
			        				ii=ii+1
#			        			w=w-1
	




#print A[a][b],a,b

	return A

###########################################################################################################################



###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções
# aqui é a geometria com base triangular


def triangular (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


	for k in range(zini,zfin+1):
		for j in range(yini,yfin+1):
		        for i in range(xini,xfin+1):

		        	a=int(i)
		        	b=int(j+(k-1)*Dy1)
		        	ii=1
		        	for h in range(1, Ntam+1):
	
		        		Npiyw=Npiy-(h-1)
		        		w=Npiyw

		        		Npixg=Npix-(h-1)
		        		g=Npixg
		        		box=0


		        		Npixgg=Npixg

		        		while (w > 0):

		        			box=box+1
		        			g=Npixgg


		        			while (g>0):

		        				div=float(radiusc/dl)




#	xa(jj)= (((2*(0)+1+(Npix-1))*radiusc/dl-i+2))
#	ya(jj)= (((2*(s-1)+1+(r-1))*radiusc/dl-j+2))
#	za(jj)= (((Sqrt(3.0)*(r-1)+1)*radiusc/dl-k+2))





		        				x[ii]= float(((float((2*(g-1)+(box)+(h-1)))*div-i+2)))
		        				y[ii]= float(((float((2*(2-2)+1+(w-1)*sqrt(3.0)+(h-1)))*div-j)+2))
		        				z[ii]= float(((float((sqrt(3.0)*(h-1)+1))*div-k+2)))





			        			if (shape=='sphere'):

				    	        		[index,aa] = incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

				    	        		if (aa==1):
					    	        		A[i][j+(k-1)*Dy1] = index



			        			if (shape=='cubic'):

				    	        		[index,aa] = incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

				    	        		if (aa==1):
					    	        		A[i][j+(k-1)*Dy1] = index



			        			if (shape=='octahedron'):	

				    	        		index = incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

				    	        		if (aa==1):
					    	        		A[i][j+(k-1)*Dy1] = index


		
		        				g=g-1


		        				ii=ii+1
		        			w=w-1
		        			Npixgg=Npixgg-1





#print A[a][b],a,b

	return A

###########################################################################################################################




###########################################################################################################################

# função ou classe que descreve qual geometria vai ser escolhida.
# primeiramente será feito com funções
# aqui é a geometria com base triangular


def triangularfcc2in2 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,cs,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


	for k in range(zini,zfin+1):
		for j in range(yini,yfin+1):
		        for i in range(xini,xfin+1):
		        	a=int(i)
		        	b=int(j+(k-1)*Dy1)
		        	ii=1		        	
		        	for o in range (1,int(Ntam/2+1)):
		        		for h in range(1,2+1):

		        			if (h==2):

#f=2*h
	
			        			f=h+1


		        				Npiyw=Npiy-(f-1)
		        				w=Npiyw
	
		        				Npixg=Npix-(f-1)
		        				g=Npixg
		        				box=0

		        				Npixgg=Npixg




		        			if(h==1):

		        				Npiyw=Npiy-(h-1)
		        				w=Npiyw

		        				Npixg=Npix-(h-1)
		        				g=Npixg
		        				box=0


		        				Npixgg=Npixg





		        			while (w > 0):

		        				box=box+1

		        				g=Npixgg



		        				while (g>0):





		        					if(h==1):


		        						x[ii]= (((2*(g-1)+(box)+(h-1))*radiusc/dl-i+2))
		        						y[ii]= (((2*(2-1)+(w-1)*sqrt(3.0)+(h-1))*radiusc/dl-j)+2)
		        						z[ii]= (((sqrt(3.0)*(h-1)+1+2*(o-1)*sqrt(3.0))*radiusc/dl-k+2))




		        					if (h==2 ):
									

#		        						print ('aqui')
		
		        						x[ii]= (((2*(g-1)+(box)+(f-1))*radiusc/dl-i+2))
		        						y[ii]= (((2*(2-1)+(w-1)*sqrt(3.0)+(f-1))*radiusc/dl-j)+2)
		        						z[ii]= (((sqrt(3.0)*(h-1)+1+2*(o-1)*sqrt(3.0))*radiusc/dl-k+2))



		        					if (shape=='sphere'):

		    				        		[index,aa] = incoreshell_sphere(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

		    				        		if (aa==1):
						    	        		A[i][j+(k-1)*Dy1] = index



		        					if (shape=='cubic'):

		    				        		[index,aa] = incoreshell_cubic(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

		    				        		if (aa==1):
						    	        		A[i][j+(k-1)*Dy1] = index



		        					if (shape=='octahedron'):	

		    				        		[index,aa] = incoreshell_octahedron(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z)

		    				        		if (aa==1):
						    	        		A[i][j+(k-1)*Dy1] = index

	

		        					if (shape=='triplate'):	

		    				        		[index,aa] = incoreshell_triplate(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber, shellnumber,dl,h,x,y,z,base,width,height,base2,width2,cs)

		    				        		if (aa==1):
						    	        		A[i][j+(k-1)*Dy1] = index



		        					if (shape=='ellipsoidal'):	

		    				        		[index,aa] = incoreshell_ellipsoid(Npiy,Npix,w,g,i,j,k,ii,a,b,cnum,radiusc,radiuscore, radiusshell, corenumber,shellnumber,dl,h,x,y,z,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)

		    				        		if (aa==1):
						    	        		A[i][j+(k-1)*Dy1] = index






		        					g=g-1





		        					ii=ii+1
		        				w=w-1
		        				Npixgg=Npixgg-1



	return A

###########################################################################################################################

def hex1 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):
#            print corenumber,shellnumber,cnum

	            a=int(i)
	            b=int(j+(k-1)*Dy1)


	            xx1 = float(((float(3.0*radiusc))/float(dl)-float(i)))  
	            yy1 = float(((float(3.0*radiusc))/float(dl)-float(j))) 
	            zz1 = float(((float(radiusc))/float(dl)-float(k+1)+float(1)))

	            xx2 = float(((float(2.0*radiusc))/float(dl)-float(i)))  
	            yy2 = float(((float((3.0-sqrt(3.0))*radiusc))/float(dl)-float(j))) 
	            zz2 = float(((float(radiusc))/float(dl)-float(k+1)+float(1)))

	            xx3 = float(((float(4.0*radiusc))/float(dl)-float(i)))  
	            yy3= float(((float((3.0-sqrt(3.0))*radiusc))/float(dl)-float(j))) 
	            zz3 = float(((float(radiusc))/float(dl)-float(k+1)+float(1)))

	            xx4 = float(((float(2.0*radiusc))/float(dl)-float(i)))  
	            yy4 = float(((float((3.0+sqrt(3.0))*radiusc))/float(dl)-float(j))) 
	            zz4 = float(((float(radiusc))/float(dl)-float(k+1)+float(1)))

	            xx5 = float(((float(4.0*radiusc))/float(dl)-float(i)))  
	            yy5 = float(((float((3.0+sqrt(3.0))*radiusc))/float(dl)-float(j))) 
	            zz5 = float(((float(radiusc))/float(dl)-float(k+1)+float(1)))

	            xx6 = float(((float(1.0*radiusc))/float(dl)-float(i)))  
	            yy6 = float(((float(3.0*radiusc))/float(dl)-float(j))) 
	            zz6 = float(((float(radiusc))/float(dl)-float(k+1)+float(1)))

	            xx7 = float(((float(5.0*radiusc))/float(dl)-float(i)))  
	            yy7 = float(((float(3.0*radiusc))/float(dl)-float(j))) 
	            zz7 = float(((float(radiusc))/float(dl)-float(k+1)+float(1)))

## second layer

	            xx8 = float(((float(4.0*radiusc))/float(dl)-float(i)))  
	            yy8 = float(((float(4.0*radiusc))/float(dl)-float(j))) 
	            zz8 = float(((float((1.0+sqrt(3.0))*radiusc))/float(dl)-float(k+1)+float(1)))

	            xx9 = float(((float(4.0*radiusc))/float(dl)-float(i)))  
	            yy9 = float(((float(2.0*radiusc))/float(dl)-float(j))) 
	            zz9 = float(((float((1.0+sqrt(3.0))*radiusc))/float(dl)-float(k+1)+float(1)))

	            xx10 = float(((float((3.0-sqrt(3.0)/2)*radiusc))/float(dl)-float(i)))  
	            yy10 = float(((float(3.0*radiusc))/float(dl)-float(j))) 
	            zz10 = float(((float((1.0+sqrt(3.0))*radiusc))/float(dl)-float(k+1)+float(1)))
	
### third layer

	            xx11 = float(((float(3.0*radiusc))/float(dl)-float(i)))  
	            yy11 = float(((float(3.0*radiusc))/float(dl)-float(j))) 
	            zz11 = float(((float((1.0+2*sqrt(3.0))*radiusc))/float(dl)-float(k+1)+float(1)))

	            a=int(i)
	            b=int(j+(k-1)*Dy1)

	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )


#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






#print A[a][b],a,b

	return A

###########################################################################################################################

###########################################################################################################################



def hex2 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):

	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):
## first layer

	            a=int(i)
	            b=int(j+(k-1)*Dy1)


        	    xx1=float((float((3)*radiusc/dl)-float(i)+1))
        	    yy1= float((float((3)*radiusc/dl)-float(j)+1))
        	    zz1=float((float(radiusc/dl)-float(k)+1))
    
        	    xx2=float((float((3-1.0)*radiusc/dl)-float(i)+1))
        	    yy2= float((float((3-sqrt(3.0)))*radiusc/dl-float(j)+1))
        	    zz2=float((float(radiusc/dl)-float(k)+1))

        	    xx3=float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
        	    yy3= float(((3.0-sqrt(3.0))*float(radiusc/dl)-j+1))
        	    zz3=float((float(radiusc/dl)-float(k)+1))

        	    xx4=float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
        	    yy4= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
        	    zz4=float((float(radiusc/dl)-float(k)+1))

        	    xx5=float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
        	    yy5= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
        	    zz5=float((float(radiusc/dl)-float(k)+1))

        	    xx6=float((1.0*float(radiusc/dl)-float(i)+1))
        	    yy6= float((3.0*float(radiusc/dl)-float(j)+1))
        	    zz6=float((float(radiusc/dl)-float(k)+1))

        	    xx7=float((5.0*float(radiusc/dl)-float(i)+1))
        	    yy7= float((3.0*float(radiusc/dl)-float(j)+1))
        	    zz7=float((float(radiusc/dl)-float(k)+1))
# second layer

        	    xx8=float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
        	    yy8= float(((3.0)*float(radiusc/dl)-float(j)+1))
        	    zz8=float(((1.0+sqrt(3.0))*float(radiusc/dl)-float(k)+1))

        	    xx9=float(((3+1)*float(radiusc/dl)-float(i)+1))
        	    yy9= float(((3.0)*float(radiusc/dl)-float(j)+1))
        	    zz9= float(((1.0+sqrt(3.0))*float(radiusc/dl)-float(k)+1))

        	    xx10= float(((3.0)*float(radiusc/dl)-float(i)+1))
        	    yy10= float(((3.0-1.0-sqrt(2.0)/2)*float(radiusc/dl)-float(j)+1))
        	    zz10= float(((1.0+sqrt(3.0))*float(radiusc/dl)-float(k)+1))

        	    xx11=float(((3.0)*float(radiusc/dl)-float(i)+1))
        	    yy11= float(((3.0+1.0+sqrt(2.0)/2)*float(radiusc/dl)-float(j)+1))
        	    zz11= float(((1.0+1.0*sqrt(3.0))*float(radiusc/dl)-float(k)+1))

# third layer

        	    xx12= float(((3.0)*float(radiusc/dl)-float(i)+1))
        	    yy12= float(((3.0-sqrt(3.0)/2)*float(radiusc/dl)-float(j)+1))
        	    zz12= float(((1.0+2*sqrt(3.0))*float(radiusc/dl)-float(k)+1))

        	    xx13= float(((3.0)*float(radiusc/dl)-float(i)+1))
        	    yy13= float(((3.0+sqrt(3.0)/2)*float(radiusc/dl)-float(j)+1))
        	    zz13= float(((1.0+2.0*sqrt(3.0))*float(radiusc/dl)-float(k)+1))

# fourth layer

        	    xx14=float(((3.0)*float(radiusc/dl)-float(i)+1))
        	    yy14= float(((3.0)*float(radiusc/dl)-float(j)+1))
        	    zz14=float(((1.0+3.0*sqrt(3.0))*float(radiusc/dl)-float(k)+1))


	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )


#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



	return A

###########################################################################################################################

###########################################################################################################################



def hex3 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):

	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):


	            a=int(i)
	            b=int(j+(k-1)*Dy1)


	            xx1=float(((3.0)*float(radiusc/dl)-float(i)+1))
	            yy1= float(((3.0)*float(radiusc/dl)-float(j)+1))
	            zz1=float((float(radiusc/dl)-float(k)+1))
    
	            xx2= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy2= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz2= float((float(radiusc/dl)-float(k)+1))

	            xx3=float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy3= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz3= float((float(radiusc/dl)-float(k)+1))

	            xx4= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy4= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz4= float((float(radiusc/dl)-float(k)+1))

	            xx5= float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy5= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz5= float((float(radiusc/dl)-float(k)+1))

	            xx6=float((1.0*float(radiusc/dl)-float(i)+1))
	            yy6= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz6= float((float(radiusc/dl)-float(k)+1))

	            xx7=float((5.0*float(radiusc/dl)-float(i)+1))
	            yy7= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz7= float((float(radiusc/dl)-float(k)+1))

# second layer


	            xx8= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy8= float(((3.0)*float(radiusc/dl)-float(j)+1))
	            zz8= float(((1.0+sqrt(3.0))*float(radiusc/dl)-float(k)+1))

	            xx9= float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy9= float(((3.0)*float(radiusc/dl)-float(j)+1))
	            zz9= float(((1+sqrt(3.0))*float(radiusc/dl)-float(k)+1))

	            xx10= float(((3.0)*float(radiusc/dl)-float(i)+1))
	            yy10= float(((3.0-1.0-sqrt(2.0)/2)*float(radiusc/dl)-float(j)+1))
	            zz10= float(((1.0+sqrt(3.0))*float(radiusc/dl)-float(k)+1))

	            xx11= float(((3.0)*float(radiusc/dl)-float(i)+1))
	            yy11= float(((3.0+1.0+sqrt(2.0)/2)*float(radiusc/dl)-float(j)+1))
	            zz11= float(((1.0+1.0*sqrt(3.0))*float(radiusc/dl)-float(k)+1))

# third layer

	            xx12= float(((3)*float(radiusc/dl)-float(i)+1))
	            yy12= float(((3)*float(radiusc/dl)-float(j)+1))
	            zz12= float(((1+2*sqrt(3.0))*float(radiusc/dl)-float(k)+1))




	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )


#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




	return A

###########################################################################################################################



###########################################################################################################################



def hex4 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):


	            a=int(i)
	            b=int(j+(k-1)*Dy1)



	            xx1= float(((3.0)*float(radiusc/dl)-float(i)+1))
	            yy1= float(((3.0)*float(radiusc/dl)-float(j)+1))
	            zz1= float((float(radiusc/dl)-float(k)+1))
    
	            xx2= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy2= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz2= float((float(radiusc/dl)-float(k)+1))

	            xx3= float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy3= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz3= float((float(radiusc/dl)-float(k)+1))

	            xx4= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy4= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz4= float((float(radiusc/dl)-float(k)+1))

	            xx5=float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy5= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz5= float((float(radiusc/dl)-float(k)+1))

	            xx6= float((1.0*float(radiusc/dl)-float(i)+1))
	            yy6= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz6= float((float(radiusc/dl)-float(k)+1))

	            xx7= float((5.0*float(radiusc/dl)-float(i)+1))
	            yy7= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz7= float((1.0*float(radiusc/dl)-float(k)+1))


#	            print(xx1,xx2,xx3,xx4,xx5,xx6,xx7)
# second layer


	            xx8= float(((3.0)*float(radiusc/dl)-float(i)+1))
	            yy8= float(((3.0)*float(radiusc/dl)-float(j)+1))
	            zz8= float((3.0*float(radiusc/dl)-float(k)+1))
    
	            xx9= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy9= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz9= float((3.0*float(radiusc/dl)-float(k)+1))

	            xx10=float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy10= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz10= float((3.0*float(radiusc/dl)-float(k)+1))

	            xx11= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy11= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz11= float((3.0*float(radiusc/dl)-float(k)+1))

	            xx12= float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy12= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz12= float((3.0*float(radiusc/dl)-float(k)+1))

	            xx13= float((1.0*float(radiusc/dl)-float(i)+1))
	            yy13= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz13= float((3.0*float(radiusc/dl)-float(k)+1))

	            xx14=float((5.0*float(radiusc/dl)-float(i)+1))
	            yy14= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz14= float((3.0*float(radiusc/dl)-float(k)+1))


# third layer


	            xx15=float(((3.0)*float(radiusc/dl)-float(i)+1))
	            yy15= float(((3.0)*float(radiusc/dl)-float(j)+1))
	            zz15= float((5.0*float(radiusc/dl)-float(k)+1))
    
	            xx16= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy16= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz16= float((5.0*float(radiusc/dl)-float(k)+1))

	            xx17= float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy17= float(((3.0-sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz17=float((5.0*float(radiusc/dl)-float(k)+1))

	            xx18= float(((3.0-1.0)*float(radiusc/dl)-float(i)+1))
	            yy18= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz18= float((5.0*float(radiusc/dl)-float(k)+1))

	            xx19= float(((3.0+1.0)*float(radiusc/dl)-float(i)+1))
	            yy19= float(((3.0+sqrt(3.0))*float(radiusc/dl)-float(j)+1))
	            zz19=float((5.0*float(radiusc/dl)-float(k)+1))

	            xx20= float((1.0*float(radiusc/dl)-float(i)+1))
	            yy20= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz20= float((5.0*float(radiusc/dl)-float(k)+1))

	            xx21= float((5.0*float(radiusc/dl)-float(i)+1))
	            yy21= float((3.0*float(radiusc/dl)-float(j)+1))
	            zz21= float((5.0*float(radiusc/dl)-float(k)+1))




	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )


#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx15/(div))**2+(yy15/(div))**2+(zz15/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx16/(div))**2+(yy16/(div))**2+(zz16/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx17/(div))**2+(yy17/(div))**2+(zz17/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx18/(div))**2+(yy18/(div))**2+(zz18/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx19/(div))**2+(yy19/(div))**2+(zz19/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx20/(div))**2+(yy20/(div))**2+(zz20/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx21/(div))**2+(yy21/(div))**2+(zz21/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum












	return A

###########################################################################################################################



###########################################################################################################################



def hex5 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):




	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):


	            a=int(i)
	            b=int(j+(k-1)*Dy1)



	            xx1=float(((3.0)*float(radiusc/dl)-float(i)+1))
	            yy1= float(((3.0)*float(radiusc/dl)-float(j)+1))
	            zz1=float((float(radiusc/dl)-k+1))
    
	            xx2=((3-1.0)*radiusc/dl-i+1)
	            yy2= ((3-sqrt(3.0))*radiusc/dl-j+1)
	            zz2=(radiusc/dl-k+1)

	            xx3=((3+1.0)*radiusc/dl-i+1)
	            yy3= ((3-sqrt(3.0))*radiusc/dl-j+1)
	            zz3=(radiusc/dl-k+1)

	            xx4=((3-1.0)*radiusc/dl-i+1)
	            yy4= ((3+sqrt(3.0))*radiusc/dl-j+1)
	            zz4=(radiusc/dl-k+1)

	            xx5=((3+1.0)*radiusc/dl-i+1)
	            yy5= ((3+sqrt(3.0))*radiusc/dl-j+1)
	            zz5=(radiusc/dl-k+1)

	            xx6=(1*radiusc/dl-i+1)
	            yy6= (3*radiusc/dl-j+1)
	            zz6=(radiusc/dl-k+1)

	            xx7=(5*radiusc/dl-i+1)
	            yy7= (3*radiusc/dl-j+1)
	            zz7=(1*radiusc/dl-k+1)

# second layer

	            xx8=((3)*radiusc/dl-i+1)
	            yy8= ((3)*radiusc/dl-j+1)
	            zz8=(3*radiusc/dl-k+1)
    
	            xx9=((3-1.0)*radiusc/dl-i+1)
	            yy9= ((3-sqrt(3.0))*radiusc/dl-j+1)
	            zz9=(3*radiusc/dl-k+1)

	            xx10=((3+1.0)*radiusc/dl-i+1)
	            yy10= ((3-sqrt(3.0))*radiusc/dl-j+1)
	            zz10=(3*radiusc/dl-k+1)

	            xx11=float(((3-1.0)*radiusc/dl-i+1))
	            yy11= float(((3+sqrt(3.0))*radiusc/dl-j+1))
	            zz11=float((3*radiusc/dl-k+1))

	            xx12=((3+1.0)*radiusc/dl-i+1)
	            yy12= ((3+sqrt(3.0))*radiusc/dl-j+1)
	            zz12=(3*radiusc/dl-k+1)

	            xx13=(1*radiusc/dl-i+1)
	            yy13= (3*radiusc/dl-j+1)
	            zz13=(3*radiusc/dl-k+1)

	            xx14=(5*radiusc/dl-i+1)
	            yy14= (3*radiusc/dl-j+1)
	            zz14=(3*radiusc/dl-k+1)

# third layer

	            xx15=((3-1)*radiusc/dl-i+1)
	            yy15= ((3)*radiusc/dl-j+1)
	            zz15=((3+sqrt(3.0))*radiusc/dl-k+1)

	            xx16=((3+1)*radiusc/dl-i+1)
	            yy16= ((3)*radiusc/dl-j+1)
	            zz16=((3+sqrt(3.0))*radiusc/dl-k+1)

	            xx17=((3)*radiusc/dl-i+1)
	            yy17= ((3-1-sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz17=((3+sqrt(3.0))*radiusc/dl-k+1)

	            xx18=((3)*radiusc/dl-i+1)
	            yy18= ((3+1+sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz18=((3+1*sqrt(3.0))*radiusc/dl-k+1)

# fourth layer

	            xx19=((3)*radiusc/dl-i+1)
	            yy19= ((3-sqrt(3.0)/2)*radiusc/dl-j+1)
	            zz19=((3+2*sqrt(3.0))*radiusc/dl-k+1)

	            xx20=((3)*radiusc/dl-i+1)
	            yy20= ((3+sqrt(3.0)/2)*radiusc/dl-j+1)
	            zz20=((3+2*sqrt(3.0))*radiusc/dl-k+1)

# fifth layer

	            xx21=((3)*radiusc/dl-i+1)
	            yy21= ((3)*radiusc/dl-j+1)
	            zz21=((4+3*sqrt(3.0))*radiusc/dl-k+1)


	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )

#	            print (xx11,yy11,zz11,3+sqrt(3.0),float(radiusc/dl))

#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

#        	        print('aqui')


        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx15/(div))**2+(yy15/(div))**2+(zz15/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx16/(div))**2+(yy16/(div))**2+(zz16/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx17/(div))**2+(yy17/(div))**2+(zz17/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx18/(div))**2+(yy18/(div))**2+(zz18/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx19/(div))**2+(yy19/(div))**2+(zz19/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx20/(div))**2+(yy20/(div))**2+(zz20/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx21/(div))**2+(yy21/(div))**2+(zz21/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




	return A




###########################################################################################################################





###########################################################################################################################



def hex6 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):




	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):


	            a=int(i)
	            b=int(j+(k-1)*Dy1)



	            xx1=((6)*radiusc/dl-i+1)
	            yy1= ((6)*radiusc/dl-j+1)
	            zz1=(radiusc/dl-k+1)
    
	            xx2=((6-1.0)*radiusc/dl-i+1)
	            yy2= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz2=(radiusc/dl-k+1)

	            xx3=((6+1.0)*radiusc/dl-i+1)
	            yy3= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz3=(radiusc/dl-k+1)

	            xx4=((6-1.0)*radiusc/dl-i+1)
	            yy4= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz4=(radiusc/dl-k+1)

	            xx5=((6+1.0)*radiusc/dl-i+1)
	            yy5= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz5=(radiusc/dl-k+1)

	            xx6=(4*radiusc/dl-i+1)
	            yy6= (6*radiusc/dl-j+1)
	            zz6=(radiusc/dl-k+1)

	            xx7=(8*radiusc/dl-i+1)
	            yy7= (6*radiusc/dl-j+1)
	            zz7=(1*radiusc/dl-k+1)

	            xx8=((6)*radiusc/dl-i+1)
	            yy8= ((8+sqrt(2.0))*radiusc/dl-j+1)
	            zz8=(1*radiusc/dl-k+1)
    
	            xx9=((6)*radiusc/dl-i+1)
	            yy9= ((4-sqrt(2.0))*radiusc/dl-j+1)
	            zz9=(1*radiusc/dl-k+1)

	            xx10=((10.0)*radiusc/dl-i+1)
	            yy10= ((6.0)*radiusc/dl-j+1)
	            zz10=(1*radiusc/dl-k+1)

	            xx11=((3-1.0)*radiusc/dl-i+1)
	            yy11= ((6.0)*radiusc/dl-j+1)
	            zz11=(1*radiusc/dl-k+1)

	            xx12=((3)*radiusc/dl-i+1)
	            yy12= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz12=(1*radiusc/dl-k+1)

	            xx13=(9*radiusc/dl-i+1)
	            yy13= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz13=(1*radiusc/dl-k+1)

	            xx14=((9+0.05)*radiusc/dl-i+1)
	            yy14= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz14=(1*radiusc/dl-k+1)


	            xx15=((3)*radiusc/dl-i+1)
	            yy15= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz15=((1)*radiusc/dl-k+1)

	            xx16=((3+1+0.05)*radiusc/dl-i+1)
	            yy16= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz16=((1)*radiusc/dl-k+1)

	            xx17=((8)*radiusc/dl-i+1)
	            yy17= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz17=((1)*radiusc/dl-k+1)

	            xx18=((8)*radiusc/dl-i+1)
	            yy18= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz18=((1)*radiusc/dl-k+1)


	            xx19=((4)*radiusc/dl-i+1)
	            yy19= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz19=((1)*radiusc/dl-k+1)

# second layer


	            xx20=((6)*radiusc/dl-i+1)
	            yy20= ((6)*radiusc/dl-j+1)
	            zz20=(3*radiusc/dl-k+1)
    
	            xx21=((6-1.0)*radiusc/dl-i+1)
	            yy21= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz21=(3*radiusc/dl-k+1)

	            xx22=((6+1.0)*radiusc/dl-i+1)
	            yy22= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz22=(3*radiusc/dl-k+1)

	            xx23=((6-1.0)*radiusc/dl-i+1)
	            yy23= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz23=(3*radiusc/dl-k+1)

	            xx24=((6+1.0)*radiusc/dl-i+1)
	            yy24= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz24=(3*radiusc/dl-k+1)

	            xx25=(4*radiusc/dl-i+1)
	            yy25= (6*radiusc/dl-j+1)
	            zz25=(3*radiusc/dl-k+1)

	            xx26=(8*radiusc/dl-i+1)
	            yy26= (6*radiusc/dl-j+1)
	            zz26=(3*radiusc/dl-k+1)

	            xx27=((6)*radiusc/dl-i+1)
	            yy27= ((8+sqrt(2.0))*radiusc/dl-j+1)
	            zz27=(3*radiusc/dl-k+1)
    
	            xx28=((6)*radiusc/dl-i+1)
	            yy28= ((4-sqrt(2.0))*radiusc/dl-j+1)
	            zz28=(3*radiusc/dl-k+1)

	            xx29=((10.0)*radiusc/dl-i+1)
	            yy29= ((6.0)*radiusc/dl-j+1)
	            zz29=(3*radiusc/dl-k+1)

	            xx30=((3-1.0)*radiusc/dl-i+1)
	            yy30= ((6.0)*radiusc/dl-j+1)
	            zz30=(3*radiusc/dl-k+1)

	            xx31=((3)*radiusc/dl-i+1)
	            yy31= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz31=(3*radiusc/dl-k+1)

	            xx32=(9*radiusc/dl-i+1)
	            yy32= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz32=(3*radiusc/dl-k+1)

	            xx33=((9+0.05)*radiusc/dl-i+1)
	            yy33= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz33=(3*radiusc/dl-k+1)


	            xx34=((3)*radiusc/dl-i+1)
	            yy34= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz34=((3)*radiusc/dl-k+1)

	            xx35=((3+1+0.05)*radiusc/dl-i+1)
	            yy35= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz35=((3)*radiusc/dl-k+1)

	            xx36=((8)*radiusc/dl-i+1)
	            yy36= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz36=((3)*radiusc/dl-k+1)

	            xx37=((8)*radiusc/dl-i+1)
	            yy37= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz37=((3)*radiusc/dl-k+1)


	            xx38=((4)*radiusc/dl-i+1)
	            yy38= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz38=((3)*radiusc/dl-k+1)


# third layer

	            xx39=((6)*radiusc/dl-i+1)
	            yy39= ((6)*radiusc/dl-j+1)
	            zz39=(5*radiusc/dl-k+1)
    
	            xx40=((6-1.0)*radiusc/dl-i+1)
	            yy40= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz40=(5*radiusc/dl-k+1)

	            xx41=((6+1.0)*radiusc/dl-i+1)
	            yy41= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz41=(5*radiusc/dl-k+1)

	            xx42=((6-1.0)*radiusc/dl-i+1)
	            yy42= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz42=(5*radiusc/dl-k+1)

	            xx43=((6+1.0)*radiusc/dl-i+1)
	            yy43= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz43=(5*radiusc/dl-k+1)

	            xx44=(4*radiusc/dl-i+1)
	            yy44= (6*radiusc/dl-j+1)
	            zz44=(5*radiusc/dl-k+1)

	            xx45=(8*radiusc/dl-i+1)
	            yy45= (6*radiusc/dl-j+1)
	            zz45=(5*radiusc/dl-k+1)

	            xx46=((6)*radiusc/dl-i+1)
	            yy46= ((8+sqrt(2.0))*radiusc/dl-j+1)
	            zz46=(5*radiusc/dl-k+1)
    
	            xx47=((6)*radiusc/dl-i+1)
	            yy47= ((4-sqrt(2.0))*radiusc/dl-j+1)
	            zz47=(5*radiusc/dl-k+1)

	            xx48=((10.0)*radiusc/dl-i+1)
	            yy48= ((6.0)*radiusc/dl-j+1)
	            zz48=(5*radiusc/dl-k+1)

	            xx49=((3-1.0)*radiusc/dl-i+1)
	            yy49= ((6.0)*radiusc/dl-j+1)
	            zz49=(5*radiusc/dl-k+1)

	            xx50=((3)*radiusc/dl-i+1)
	            yy50= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz50=(5*radiusc/dl-k+1)

	            xx51=(9*radiusc/dl-i+1)
	            yy51= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz51=(5*radiusc/dl-k+1)

	            xx52=((9+0.05)*radiusc/dl-i+1)
	            yy52= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz52=(5*radiusc/dl-k+1)


	            xx53=((3)*radiusc/dl-i+1)
	            yy53= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz53=((5)*radiusc/dl-k+1)

	            xx54=((3+1+0.05)*radiusc/dl-i+1)
	            yy54= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz54=((5)*radiusc/dl-k+1)

	            xx55=((8)*radiusc/dl-i+1)
	            yy55= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz55=((5)*radiusc/dl-k+1)

	            xx56=((8)*radiusc/dl-i+1)
	            yy56= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz56=((5)*radiusc/dl-k+1)


	            xx57=((4)*radiusc/dl-i+1)
	            yy57= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz57=((5)*radiusc/dl-k+1)




	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )

#	            print (xx11,yy11,zz11,3+sqrt(3.0),float(radiusc/dl))

#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

#        	        print('aqui')


        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx15/(div))**2+(yy15/(div))**2+(zz15/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx16/(div))**2+(yy16/(div))**2+(zz16/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx17/(div))**2+(yy17/(div))**2+(zz17/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx18/(div))**2+(yy18/(div))**2+(zz18/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx19/(div))**2+(yy19/(div))**2+(zz19/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx20/(div))**2+(yy20/(div))**2+(zz20/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx21/(div))**2+(yy21/(div))**2+(zz21/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx22/(div))**2+(yy22/(div))**2+(zz22/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx23/(div))**2+(yy23/(div))**2+(zz23/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx24/(div))**2+(yy24/(div))**2+(zz24/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx25/(div))**2+(yy25/(div))**2+(zz25/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx26/(div))**2+(yy26/(div))**2+(zz26/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum








        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx27/(div))**2+(yy27/(div))**2+(zz27/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx28/(div))**2+(yy28/(div))**2+(zz28/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx29/(div))**2+(yy29/(div))**2+(zz29/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx30/(div))**2+(yy30/(div))**2+(zz30/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx31/(div))**2+(yy31/(div))**2+(zz31/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx32/(div))**2+(yy32/(div))**2+(zz32/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx33/(div))**2+(yy33/(div))**2+(zz33/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx34/(div))**2+(yy34/(div))**2+(zz34/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx35/(div))**2+(yy35/(div))**2+(zz35/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx36/(div))**2+(yy36/(div))**2+(zz36/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx37/(div))**2+(yy37/(div))**2+(zz37/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx38/(div))**2+(yy38/(div))**2+(zz38/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx39/(div))**2+(yy39/(div))**2+(zz39/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx40/(div))**2+(yy40/(div))**2+(zz40/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx41/(div))**2+(yy41/(div))**2+(zz41/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx42/(div))**2+(yy42/(div))**2+(zz42/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx43/(div))**2+(yy43/(div))**2+(zz43/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx44/(div))**2+(yy44/(div))**2+(zz44/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum








        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx45/(div))**2+(yy45/(div))**2+(zz45/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx46/(div))**2+(yy46/(div))**2+(zz46/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx47/(div))**2+(yy47/(div))**2+(zz47/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx48/(div))**2+(yy48/(div))**2+(zz48/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx49/(div))**2+(yy49/(div))**2+(zz49/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx50/(div))**2+(yy50/(div))**2+(zz50/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx51/(div))**2+(yy51/(div))**2+(zz51/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx52/(div))**2+(yy52/(div))**2+(zz52/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx53/(div))**2+(yy53/(div))**2+(zz53/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx54/(div))**2+(yy54/(div))**2+(zz54/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx55/(div))**2+(yy55/(div))**2+(zz55/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx56/(div))**2+(yy56/(div))**2+(zz56/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx57/(div))**2+(yy57/(div))**2+(zz57/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




	return A




###########################################################################################################################


###########################################################################################################################






def hex7 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):


	            a=int(i)
	            b=int(j+(k-1)*Dy1)




	            xx1=((6)*radiusc/dl-i+1)
	            yy1= ((6)*radiusc/dl-j+1)
	            zz1=(radiusc/dl-k+1)
    
	            xx2=((6-1.0)*radiusc/dl-i+1)
	            yy2= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz2=(radiusc/dl-k+1)

	            xx3=((6+1.0)*radiusc/dl-i+1)
	            yy3= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz3=(radiusc/dl-k+1)

	            xx4=((6-1.0)*radiusc/dl-i+1)
	            yy4= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz4=(radiusc/dl-k+1)

	            xx5=((6+1.0)*radiusc/dl-i+1)
	            yy5= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz5=(radiusc/dl-k+1)

	            xx6=(4*radiusc/dl-i+1)
	            yy6= (6*radiusc/dl-j+1)
	            zz6=(radiusc/dl-k+1)

	            xx7=(8*radiusc/dl-i+1)
	            yy7= (6*radiusc/dl-j+1)
	            zz7=(1*radiusc/dl-k+1)

	            xx8=((6)*radiusc/dl-i+1)
	            yy8= ((8+sqrt(2.0))*radiusc/dl-j+1)
	            zz8=(1*radiusc/dl-k+1)
    

	            xx9=((6)*radiusc/dl-i+1)
	            yy9= ((4-sqrt(2.0))*radiusc/dl-j+1)
	            zz9=(1*radiusc/dl-k+1)

	            xx10=((10.0)*radiusc/dl-i+1)
	            yy10= ((6.0)*radiusc/dl-j+1)
	            zz10=(1*radiusc/dl-k+1)

	            xx11=((3-1.0)*radiusc/dl-i+1)
	            yy11= ((6.0)*radiusc/dl-j+1)
	            zz11=(1*radiusc/dl-k+1)

	            xx12=((3)*radiusc/dl-i+1)
	            yy12= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz12=(1*radiusc/dl-k+1)

	            xx13=(9*radiusc/dl-i+1)
	            yy13= ((6-sqrt(3.0))*radiusc/dl-j+1)
	            zz13=(1*radiusc/dl-k+1)

	            xx14=((9+0.05)*radiusc/dl-i+1)
	            yy14= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz14=(1*radiusc/dl-k+1)


	            xx15=((3)*radiusc/dl-i+1)
	            yy15= ((6+sqrt(3.0))*radiusc/dl-j+1)
	            zz15=((1)*radiusc/dl-k+1)


	            xx16=((3+1+0.05)*radiusc/dl-i+1)
	            yy16= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz16=((1)*radiusc/dl-k+1)


	            xx17=((8)*radiusc/dl-i+1)
	            yy17= ((11-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz17=((1)*radiusc/dl-k+1)



	            xx18=((8)*radiusc/dl-i+1)
	            yy18= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz18=((1)*radiusc/dl-k+1)


	            xx19=((4)*radiusc/dl-i+1)
	            yy19= ((1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz19=((1)*radiusc/dl-k+1)

# second layer


	            xx20=((5)*radiusc/dl-i+1)
	            yy20= ((6)*radiusc/dl-j+1)
	            zz20=((1+sqrt(3.0))*radiusc/dl-k+1)
    
	            xx21=((7.0)*radiusc/dl-i+1)
	            yy21= ((6.0)*radiusc/dl-j+1)
	            zz21=((1+sqrt(3.0))*radiusc/dl-k+1)



	            xx22=((6.0)*radiusc/dl-i+1)
	            yy22= ((7.0+sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz22=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx23=((6)*radiusc/dl-i+1)
	            yy23= ((5-sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz23=((1+sqrt(3.0))*radiusc/dl-k+1)

# third layer

	            xx24=((6.0)*radiusc/dl-i+1)
	            yy24= ((6.0-sqrt(3.0)/2)*radiusc/dl-j+1)
	            zz24=((1+2*sqrt(3.0))*radiusc/dl-k+1)


	            xx25=((6)*radiusc/dl-i+1)
	            yy25= ((6+sqrt(3.0)/2)*radiusc/dl-j+1)
	            zz25=((1+2*sqrt(3.0))*radiusc/dl-k+1)

# fourth layer

	            xx26=((6.0)*radiusc/dl-i+1)
	            yy26= ((6.0)*radiusc/dl-j+1)
	            zz26=((1+3*sqrt(3.0))*radiusc/dl-k+1)








	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )

#	            print (xx11,yy11,zz11,3+sqrt(3.0),float(radiusc/dl))

#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

#        	        print('aqui')


        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx15/(div))**2+(yy15/(div))**2+(zz15/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx16/(div))**2+(yy16/(div))**2+(zz16/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx17/(div))**2+(yy17/(div))**2+(zz17/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx18/(div))**2+(yy18/(div))**2+(zz18/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx19/(div))**2+(yy19/(div))**2+(zz19/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx20/(div))**2+(yy20/(div))**2+(zz20/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx21/(div))**2+(yy21/(div))**2+(zz21/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx22/(div))**2+(yy22/(div))**2+(zz22/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx23/(div))**2+(yy23/(div))**2+(zz23/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx24/(div))**2+(yy24/(div))**2+(zz24/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx25/(div))**2+(yy25/(div))**2+(zz25/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx26/(div))**2+(yy26/(div))**2+(zz26/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






	return A

##################################################################################################################################



###########################################################################################################################






def hex8 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):


	            a=int(i)
	            b=int(j+(k-1)*Dy1)






	            xx1=((6)*radiusc/dl-i+1)
	            yy1= ((6+1)*radiusc/dl-j+1)
	            zz1=(radiusc/dl-k+1)
    
	            xx2=((6-1.0)*radiusc/dl-i+1)
	            yy2= ((6+1-sqrt(3.0))*radiusc/dl-j+1)
	            zz2=(radiusc/dl-k+1)

	            xx3=((6+1.0)*radiusc/dl-i+1)
	            yy3= ((6+1-sqrt(3.0))*radiusc/dl-j+1)
	            zz3=(radiusc/dl-k+1)

	            xx4=((6-1.0)*radiusc/dl-i+1)
	            yy4= ((6+1+sqrt(3.0))*radiusc/dl-j+1)
	            zz4=(radiusc/dl-k+1)

	            xx5=((6+1.0)*radiusc/dl-i+1)
	            yy5= ((6+1+sqrt(3.0))*radiusc/dl-j+1)
	            zz5=(radiusc/dl-k+1)

	            xx6=(4*radiusc/dl-i+1)
	            yy6= ((6+1)*radiusc/dl-j+1)
	            zz6=(radiusc/dl-k+1)

	            xx7=(8*radiusc/dl-i+1)
	            yy7= ((6+1)*radiusc/dl-j+1)
	            zz7=(1*radiusc/dl-k+1)

	            xx8=((6)*radiusc/dl-i+1)
	            yy8= ((8+1+sqrt(2.0))*radiusc/dl-j+1)
	            zz8=(1*radiusc/dl-k+1)
    
	            xx9=((6)*radiusc/dl-i+1)
	            yy9= ((4+1-sqrt(2.0))*radiusc/dl-j+1)
	            zz9=(1*radiusc/dl-k+1)

	            xx10=((10.0)*radiusc/dl-i+1)
	            yy10= ((6.0+1)*radiusc/dl-j+1)
	            zz10=(1*radiusc/dl-k+1)

	            xx11=((3-1.0)*radiusc/dl-i+1)
	            yy11= ((6.0+1)*radiusc/dl-j+1)
	            zz11=(1*radiusc/dl-k+1)

	            xx12=((3)*radiusc/dl-i+1)
	            yy12= ((6+1-sqrt(3.0))*radiusc/dl-j+1)
	            zz12=(1*radiusc/dl-k+1)

	            xx13=(9*radiusc/dl-i+1)
	            yy13= ((6+1-sqrt(3.0))*radiusc/dl-j+1)
	            zz13=(1*radiusc/dl-k+1)

	            xx14=((9+0.05)*radiusc/dl-i+1)
	            yy14= ((6+1+sqrt(3.0))*radiusc/dl-j+1)
	            zz14=(1*radiusc/dl-k+1)


	            xx15=((3)*radiusc/dl-i+1)
	            yy15= ((6+1+sqrt(3.0))*radiusc/dl-j+1)
	            zz15=((1)*radiusc/dl-k+1)

	            xx16=((3+1+0.05)*radiusc/dl-i+1)
	            yy16= ((11+1-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz16=((1)*radiusc/dl-k+1)

	            xx17=((8)*radiusc/dl-i+1)
	            yy17= ((11+1-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz17=((1)*radiusc/dl-k+1)

	            xx18=((8)*radiusc/dl-i+1)
	            yy18= ((1+1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz18=((1)*radiusc/dl-k+1)


	            xx19=((4)*radiusc/dl-i+1)
	            yy19= ((1+1+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz19=((1)*radiusc/dl-k+1)

# second layer


	            xx20 = (6*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy20 = (7*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz20 = (3*(radiuscore+(radiusc-radiuscore))/dl-k+1);




	            xx21 = (5*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy21 = ((7-sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz21 = (3*(radiuscore+(radiusc-radiuscore))/dl-k+1);


	            xx22 = (7*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy22 = ((7-sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz22 = (3*(radiuscore+(radiusc-radiuscore))/dl-k+1);




	            xx23 = (5*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy23 = ((7+sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz23 = (3*(radiuscore+(radiusc-radiuscore))/dl-k+1);




	            xx24 = (7*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy24 = ((7+sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz24 = (3*(radiuscore+(radiusc-radiuscore))/dl-k+1);



	            xx25 = (4*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy25 = (7*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz25 = (3*(radiuscore+(radiusc-radiuscore))/dl-k+1);




	            xx26 = (8*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy26 = (7*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz26 = (3*(radiuscore+(radiusc-radiuscore))/dl-k+1);


# third layer

	            xx27 = (5*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy27 = (7*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz27 = ((3+sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-k+1);


	            xx28 = (7*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy28 = (7*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz28 = ((3+sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-k+1);


	            xx29 = (6*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy29 = ((8.0+sqrt(2.0)/2)*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz29 = ((3+sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-k+1);




	            xx30 = (6*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy30 = ((6.0-sqrt(2.0)/2)*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz30 = ((3+sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-k+1);


# fourth layer

	            xx31 = (6*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy31 = ((7.0-sqrt(3.0)/2)*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz31 = ((3+2*sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-k+1);


	            xx32 = (6*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy32 = ((7.0+sqrt(3.0)/2)*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz32 = ((3+2*sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-k+1);

# fifth layer


	            xx33 = (6*(radiuscore + (radiusc-radiuscore))/dl-i); 
	            yy33 = ((7.0)*(radiuscore+(radiusc-radiuscore))/dl-j); 
	            zz33 = ((3+3*sqrt(3.0))*(radiuscore+(radiusc-radiuscore))/dl-k+1);








	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )

#	            print (xx11,yy11,zz11,3+sqrt(3.0),float(radiusc/dl))

#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

#        	        print('aqui')


        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx15/(div))**2+(yy15/(div))**2+(zz15/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx16/(div))**2+(yy16/(div))**2+(zz16/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx17/(div))**2+(yy17/(div))**2+(zz17/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx18/(div))**2+(yy18/(div))**2+(zz18/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx19/(div))**2+(yy19/(div))**2+(zz19/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx20/(div))**2+(yy20/(div))**2+(zz20/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx21/(div))**2+(yy21/(div))**2+(zz21/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx22/(div))**2+(yy22/(div))**2+(zz22/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx23/(div))**2+(yy23/(div))**2+(zz23/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx24/(div))**2+(yy24/(div))**2+(zz24/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx25/(div))**2+(yy25/(div))**2+(zz25/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx26/(div))**2+(yy26/(div))**2+(zz26/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum








        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx27/(div))**2+(yy27/(div))**2+(zz27/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx28/(div))**2+(yy28/(div))**2+(zz28/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx29/(div))**2+(yy29/(div))**2+(zz29/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx30/(div))**2+(yy30/(div))**2+(zz30/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx31/(div))**2+(yy31/(div))**2+(zz31/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx32/(div))**2+(yy32/(div))**2+(zz32/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx33/(div))**2+(yy33/(div))**2+(zz33/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







	return A

##################################################################################################################################


###########################################################################################################################






def hex9 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs):


	for k in range(zini,zfin+1):
	    for j in range(yini,yfin+1):
	        for i in range(xini,xfin+1):


	            a=int(i)
	            b=int(j+(k-1)*Dy1)



	            xx1=((6)*radiusc/dl-i+1)
	            yy1= ((7)*radiusc/dl-j+1)
	            zz1=(radiusc/dl-k+1)
    
	            xx2=((6-1.0)*radiusc/dl-i+1)
	            yy2= ((7-sqrt(3.0))*radiusc/dl-j+1)
	            zz2=(radiusc/dl-k+1)

	            xx3=((6+1.0)*radiusc/dl-i+1)
	            yy3= ((7-sqrt(3.0))*radiusc/dl-j+1)
	            zz3=(radiusc/dl-k+1)

	            xx4=((6-1.0)*radiusc/dl-i+1)
	            yy4= ((7+sqrt(3.0))*radiusc/dl-j+1)
	            zz4=(radiusc/dl-k+1)

	            xx5=((6+1.0)*radiusc/dl-i+1)
	            yy5= ((7+sqrt(3.0))*radiusc/dl-j+1)
	            zz5=(radiusc/dl-k+1)

	            xx6=(4*radiusc/dl-i+1)
	            yy6= (7*radiusc/dl-j+1)
	            zz6=(radiusc/dl-k+1)

	            xx7=(8*radiusc/dl-i+1)
	            yy7= (7*radiusc/dl-j+1)
	            zz7=(1*radiusc/dl-k+1)

	            xx8=((6)*radiusc/dl-i+1)
	            yy8= ((9+sqrt(2.0))*radiusc/dl-j+1)
	            zz8=(1*radiusc/dl-k+1)
    
	            xx9=((6)*radiusc/dl-i+1)
	            yy9= ((5-sqrt(2.0))*radiusc/dl-j+1)
	            zz9=(1*radiusc/dl-k+1)

	            xx10=((10.0)*radiusc/dl-i+1)
	            yy10= ((7.0)*radiusc/dl-j+1)
	            zz10=(1*radiusc/dl-k+1)

	            xx11=((3-1.0)*radiusc/dl-i+1)
	            yy11= ((7.0)*radiusc/dl-j+1)
	            zz11=(1*radiusc/dl-k+1)

	            xx12=((3)*radiusc/dl-i+1)
	            yy12= ((7-sqrt(3.0))*radiusc/dl-j+1)
	            zz12=(1*radiusc/dl-k+1)

	            xx13=(9*radiusc/dl-i+1)
	            yy13= ((7-sqrt(3.0))*radiusc/dl-j+1)
	            zz13=(1*radiusc/dl-k+1)

	            xx14=((9+0.05)*radiusc/dl-i+1)
	            yy14= ((7+sqrt(3.0))*radiusc/dl-j+1)
	            zz14=(1*radiusc/dl-k+1)


	            xx15=((3)*radiusc/dl-i+1)
	            yy15= ((7+sqrt(3.0))*radiusc/dl-j+1)
	            zz15=((1)*radiusc/dl-k+1)

	            xx16=((3+1+0.05)*radiusc/dl-i+1)
	            yy16= ((12-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz16=((1)*radiusc/dl-k+1)

	            xx17=((8)*radiusc/dl-i+1)
	            yy17= ((12-sqrt(3.0)+0.15)*radiusc/dl-j+1)
	            zz17=((1)*radiusc/dl-k+1)

	            xx18=((8)*radiusc/dl-i+1)
	            yy18= ((2+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz18=((1)*radiusc/dl-k+1)


	            xx19=((4)*radiusc/dl-i+1)
	            yy19= ((2+sqrt(3.0)-0.1)*radiusc/dl-j+1)
	            zz19=((1)*radiusc/dl-k+1)

# second layer

	            xx20=((5)*radiusc/dl-i+1)
	            yy20= ((7)*radiusc/dl-j+1)
	            zz20=((1+sqrt(3.0))*radiusc/dl-k+1)
    
	            xx21=((7)*radiusc/dl-i+1)
	            yy21= ((7)*radiusc/dl-j+1)
	            zz21=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx22=((6)*radiusc/dl-i+1)
	            yy22= ((8+sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz22=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx23=((6)*radiusc/dl-i+1)
	            yy23= ((6-sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz23=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx24=((3.0)*radiusc/dl-i+1)
	            yy24= ((7)*radiusc/dl-j+1)
	            zz24=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx25=(9*radiusc/dl-i+1)
	            yy25= (7*radiusc/dl-j+1)
	            zz25=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx26=((8)*radiusc/dl-i+1)
	            yy26= ((8+sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz26=((1+sqrt(3.0))*radiusc/dl-k+1)



	            xx27=((4)*radiusc/dl-i+1)
	            yy27= ((8+sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz27=((1+sqrt(3.0))*radiusc/dl-k+1)
    
	            xx28=((8.0)*radiusc/dl-i+1)
	            yy28= ((6.0-sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz28=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx29=((4.0)*radiusc/dl-i+1)
	            yy29= ((6.0-sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz29=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx30=((5)*radiusc/dl-i+1)
	            yy30= ((8+sqrt(2.0)/2+sqrt(3.0))*radiusc/dl-j+1)
	            zz30=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx31=((7.0)*radiusc/dl-i+1)
	            yy31= ((8.0+sqrt(2.0)/2+sqrt(3.0))*radiusc/dl-j+1)
	            zz31=((1+sqrt(3.0))*radiusc/dl-k+1)

	            xx32=((5)*radiusc/dl-i+1)
	            yy32= ((6-sqrt(2.0)/2-sqrt(3.0))*radiusc/dl-j+1)
	            zz32=((1+sqrt(3.0))*radiusc/dl-k+1)


	            xx33=((7.0)*radiusc/dl-i+1)
	            yy33= ((6.0-sqrt(2.0)/2-sqrt(3.0))*radiusc/dl-j+1)
	            zz33=((1+sqrt(3.0))*radiusc/dl-k+1)

# third layer


	            xx34=((4.0)*radiusc/dl-i+1)
	            yy34= ((7.0)*radiusc/dl-j+1)
	            zz34=((1+2*sqrt(3.0))*radiusc/dl-k+1)

	            xx35=((6)*radiusc/dl-i+1)
	            yy35= ((7)*radiusc/dl-j+1)
	            zz35=((1+2*sqrt(3.0))*radiusc/dl-k+1)


	            xx36=((8.0)*radiusc/dl-i+1)
	            yy36= ((7.0)*radiusc/dl-j+1)
	            zz36=((1+2*sqrt(3.0))*radiusc/dl-k+1)


	            xx37=((7.0)*radiusc/dl-i+1)
	            yy37= ((8.0+sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz37=((1+2*sqrt(3.0))*radiusc/dl-k+1)

	            xx38=((5)*radiusc/dl-i+1)
	            yy38= ((8+sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz38=((1+2*sqrt(3.0))*radiusc/dl-k+1)


	            xx39=((7.0)*radiusc/dl-i+1)
	            yy39= ((6.0-sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz39=((1+2*sqrt(3.0))*radiusc/dl-k+1)

	            xx40=((5.0)*radiusc/dl-i+1)
	            yy40= ((6.0-sqrt(2.0)/2)*radiusc/dl-j+1)
	            zz40=((1+2*sqrt(3.0))*radiusc/dl-k+1)

	            xx41=((6)*radiusc/dl-i+1)
	            yy41= ((8+sqrt(2.0)/2+sqrt(3.0))*radiusc/dl-j+1)
	            zz41=((1+2*sqrt(3.0))*radiusc/dl-k+1)


	            xx42=((6.0)*radiusc/dl-i+1)
	            yy42= ((6.0-sqrt(2.0)/2-sqrt(3.0))*radiusc/dl-j+1)
	            zz42=((1+2*sqrt(3.0))*radiusc/dl-k+1)

# fourth layer


	            xx43=((5)*radiusc/dl-i+1)
	            yy43= ((7.0)*radiusc/dl-j+1)
	            zz43=((1+3*sqrt(3.0))*radiusc/dl-k+1)



	            xx44=((7.0)*radiusc/dl-i+1)
	            yy44= ((7.0)*radiusc/dl-j+1)
	            zz44=((1+3*sqrt(3.0))*radiusc/dl-k+1)

	            xx45=((6)*radiusc/dl-i+1)
	            yy45= ((7+sqrt(3.0))*radiusc/dl-j+1)
	            zz45=((1+3*sqrt(3.0))*radiusc/dl-k+1)


	            xx46=((6.0)*radiusc/dl-i+1)
	            yy46= ((7-sqrt(3.0))*radiusc/dl-j+1)
	            zz46=((1+3*sqrt(3.0))*radiusc/dl-k+1)

# fifth layer

	            xx47=((6.0)*radiusc/dl-i+1)
	            yy47= ((7.0)*radiusc/dl-j+1)
	            zz47=((1+4*sqrt(3.0))*radiusc/dl-k+1)





	            div=float(radiusc)/float(dl)
	            cc=float(((xx1/(div))**2+(yy1/(div))**2+(zz1/(div))**2) )

#	            print (xx11,yy11,zz11,3+sqrt(3.0),float(radiusc/dl))

#            print radiuscore,radiusshell,radiusc
            
	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
	                A[a][b]=shellnumber

	            if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
	                A[a][b]=cnum

	            div=float(radiusc)/float(dl)
	            cc=float(((xx2/(div))**2+(yy2/(div))**2+(zz2/(div))**2) )


	            if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
	                A[a][b]=corenumber

	            if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx3/(div))**2+(yy3/(div))**2+(zz3/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx4/(div))**2+(yy4/(div))**2+(zz4/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx5/(div))**2+(yy5/(div))**2+(zz5/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx6/(div))**2+(yy6/(div))**2+(zz6/(div))**2) )

        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum
	

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx7/(div))**2+(yy7/(div))**2+(zz7/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx8/(div))**2+(yy8/(div))**2+(zz8/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx9/(div))**2+(yy9/(div))**2+(zz9/(div))**2) )

	
        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx10/(div))**2+(yy10/(div))**2+(zz10/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx11/(div))**2+(yy11/(div))**2+(zz11/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

#        	        print('aqui')


        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx12/(div))**2+(yy12/(div))**2+(zz12/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx13/(div))**2+(yy13/(div))**2+(zz13/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx14/(div))**2+(yy14/(div))**2+(zz14/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum



        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx15/(div))**2+(yy15/(div))**2+(zz15/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx16/(div))**2+(yy16/(div))**2+(zz16/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx17/(div))**2+(yy17/(div))**2+(zz17/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx18/(div))**2+(yy18/(div))**2+(zz18/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx19/(div))**2+(yy19/(div))**2+(zz19/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx20/(div))**2+(yy20/(div))**2+(zz20/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx21/(div))**2+(yy21/(div))**2+(zz21/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx22/(div))**2+(yy22/(div))**2+(zz22/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum

        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx23/(div))**2+(yy23/(div))**2+(zz23/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx24/(div))**2+(yy24/(div))**2+(zz24/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx25/(div))**2+(yy25/(div))**2+(zz25/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx26/(div))**2+(yy26/(div))**2+(zz26/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum








        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx27/(div))**2+(yy27/(div))**2+(zz27/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx28/(div))**2+(yy28/(div))**2+(zz28/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx29/(div))**2+(yy29/(div))**2+(zz29/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx30/(div))**2+(yy30/(div))**2+(zz30/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx31/(div))**2+(yy31/(div))**2+(zz31/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx32/(div))**2+(yy32/(div))**2+(zz32/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx33/(div))**2+(yy33/(div))**2+(zz33/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx34/(div))**2+(yy34/(div))**2+(zz34/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx35/(div))**2+(yy35/(div))**2+(zz35/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx36/(div))**2+(yy36/(div))**2+(zz36/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx37/(div))**2+(yy37/(div))**2+(zz37/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx38/(div))**2+(yy38/(div))**2+(zz38/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum


        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx39/(div))**2+(yy39/(div))**2+(zz39/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx40/(div))**2+(yy40/(div))**2+(zz40/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx41/(div))**2+(yy41/(div))**2+(zz41/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx42/(div))**2+(yy42/(div))**2+(zz42/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum




        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx43/(div))**2+(yy43/(div))**2+(zz43/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum







        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx44/(div))**2+(yy44/(div))**2+(zz44/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum








        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx45/(div))**2+(yy45/(div))**2+(zz45/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx46/(div))**2+(yy46/(div))**2+(zz46/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum






        	    div=float(radiusc)/float(dl)
        	    cc=float(((xx47/(div))**2+(yy47/(div))**2+(zz47/(div))**2) )


        	    if ( ((cc))<= (((radiuscore)/(radiuscore+(radiusc-radiuscore)) ))**2) :
        	        A[a][b]=corenumber

        	    if ( cc >= (((radiuscore)/(radiuscore+(radiusc-radiuscore)))**2) and (cc) <= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2)):
        	        A[a][b]=shellnumber

        	    if ( cc >= (((radiusshell)/(radiuscore+(radiusc-radiuscore)))**2) and cc <= 1.0):
        	        A[a][b]=cnum





	return A

##################################################################################################################################




### função que descreve a compactação


def compac(A,Dx1,Dy1,Dz1,zini,zfin,yini,yfin,xini,xfin):

#	print (A)

	print(A[21][14],"aqui agora")

	fileemp=open('/home/vagner/Desktop/GREAT_REFS/MEISMATHEMATICA/shay/matrizpython/marmittexamples/ajudamarmitt/parallelmtxgen/emp.mtx','w')

	fileemp.write(str(Dx1)+ ' ')
	fileemp.write(str(Dy1)+ ' ')
	fileemp.write(str(Dz1)+ '\n')
			
	print(A[1][1],A[2][1],A[3][1],A[4][1],A[5][1],A[6][1],A[7][1],A[21][14],'first matrix element')
			
	iCp_old = A[0][0]
	c = 0
	for k in range(zini,zfin+1):
		for j in range(yini,yfin+1):
			for i in range(xini,xfin+1):
				iCp = A[i][j+(k-1)*Dy1]
				if iCp == iCp_old:
					c += 1


#					print(A[21][14],"aqui agora porra")


					if c == 50000:
						print(c, iCp_old, file=fileemp)
						c = 1
						iCp_old = A[i][j+(k-1)*Dy1]



				else:

					print(c, iCp_old, file=fileemp)
					iCp_old = iCp
					c = 1
#					print("veio ate aqui!!!!")					
				if (i==Dx1 and j==Dy1 and k==Dz1):
				
					print(c,iCp_old, file=fileemp)
				




###########################################################################################################################


# Descreve as dimensões da matriz a ser criada, juntamente com seus parâmetros
# parâmetros que serão dados pelo usuário


#if ( Npix > Npiy):

#    Nmaior=Npix
#    Nmenor=Npiy
#elif(Npix <= Npiy):

#    Nmaior=Npiy
#    Nmenor=Npix




#	print(1)


def mtxgen(Npix, Npiy,Ntam,dl,cs,radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum,escolha,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,rcarbonpt,radiuspt,rcarbonpd,radiuspd,distcentro,alloynumber):


#	Npix=3
#	Npiy=3
#	Ntam=3

#	dl=5.0


#alloytudo= (input('Alloy de PtPd em tudo? ("s" OR "n") '))


#cascaptpd= (input('Casca de PtPd ("s" OR "n") '))



#core= (input('PT-core OR PD-core("pt" OR "pd") '))

#cs=10.0



#if(alloytudo =='s'):
#	[radiusptpd,radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum]= estruturas(alloytudo,cascaptpd,core,cs,dl)

#if(alloytudo == 'n' and cascaptpd =='n' and (core =='pd' or core == 'pt')):
#	[radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum] = estruturas(alloytudo,cascaptpd,core,cs,dl)
#[radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum]

	if(escolha=='hexmono' or escolha== 'hexpy'):
		Npix=Npix*2
		Npiy=Npiy*2-1
		print('sim')



#	if(escolha != 'triplate'):



	if(escolha != 'interpene'):

		print(Npix,Npiy)
		dx=Npix*2*radiusc
		dy=Npix*2*radiusc
		dz=Ntam*2*radiusc

		Dx1=int(dx/dl)
		Dy1=int(dy/dl)
		Dz1=int(dz/dl)

	elif(escolha=='interpene'):


		dx=radiuspt+radiuspd+(rcarbonpt-radiuspt)+(rcarbonpd-radiuspd)+10+distcentro

		if (radiuspt>=radiuspd):

			dy=2*(rcarbonpt)+distcentro
#			dz=2*radiuspt+2*(rcarbonpt-radiuspt)
			dz=dy
			dx=dy
			radiusc= rcarbonpt+distcento

		if (radiuspd>radiuspt):

			dy=2*(rcarbonpd)+distcentro
#			dz=2*radiuspd+2*(rcarbonpd-radiuspd)
			dz=dy
			dx=dy
			radiusc= rcarbonpt+distcento
	


#	elif (escolha=='triplate')



### chamar a função que zera a matriz


#zeromatrix(Dy1,Dy1,Dz1,A)
#!!!
	A = zeromatrix(Dx1, Dy1, Dz1)
#!!!

#	print(1)


	Npixn=Dx1
	Npiyn=Dy1
	Ntamn=Dz1


	xini=1
	yini=1
	zini=1
	xfin=int(Npixn)#+1)
	yfin=int(Npiyn)#+1)
	zfin=int(Ntamn)#+1)




## zera-se o array que será utilizado

#zeroarray(Nmaior,Nmenor)
#!!!
	[x, y, z, Npart] = zeroarray(Npix,Npiy,Ntam)
# escolha da geomtria
#!!!

#	escolha = 'py'


	if (escolha == 'jonder1'):
	


		A=jonderstruct1 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs)


	if (escolha == 'jonder2'):

		A=jonderstruct2 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs)



	if (escolha == 'tri'):
	
#	if(alloytudo =='s'):
#	[radiusptpd,radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum]= estruturas(alloytudo,cascaptpd,core,cs,dl)

		A=triangular (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)

	if (escolha == 'rectangular'):


#		if __name__ == '__main__':
#			pool = Pool(processes=4)              # process per core
#			pool.map(rectangular,)  # proces data_inp


#		return print(1)	
#		processes=[]
#		p = Process(target=rectangular, args=(xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs))
#		p.start()
#		processes.append(p)
		
#		print('aqui')
#		for p in processes:
#			p.join()




		A= rectangular(xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs)





	if (escolha == 'column'):

		A=columnpiling (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs)



	if (escolha == 'pygen'):
	
#	if(alloytudo =='s'):
#	[radiusptpd,radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum]= estruturas(alloytudo,cascaptpd,core,cs,dl)
#		processes=[]
#		p = Process(target=pygeneral, args=((xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs,rcarbonpt,radiuspt,rcarbonpd,radiuspd,distcentro,alloynumber)))
#		p.start()
#		processes.append(p)
		
#		print('aqui')
		for p in processes:
			p.join()

	
		A=pygeneral (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs,cs,rcarbonpt,radiuspt,rcarbonpd,radiuspd,distcentro,alloynumber)
#	!!!

	if (escolha == 'hexmono'):
	
#	if(alloytudo =='s'):
#	[radiusptpd,radiusc,radiuscore,radiusshell,corenumber,shellnumber,cnum]= estruturas(alloytudo,cascaptpd,core,cs,dl)

	
		A=hexmono (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)
#	!!!

	if (escolha == 'hexpy'):


		A=hexpyramid (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)
	
	if (escolha == 'hex1'):

		A=hex1 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)


	if (escolha == 'hex2'):

		A=hex2 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)


	if (escolha == 'hex3'):

		A=hex3 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)


	if (escolha == 'hex4'):


		A=hex4 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)



	if (escolha == 'hex5'):


		A=hex5 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)




	if (escolha == 'hex6'):


		A=hex6 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)


	if (escolha == 'hex7'):


		A=hex7 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)




	if (escolha == 'hex8'):


		A=hex8 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)


	if (escolha == 'hex9'):


		A=hex9 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)



	if (escolha == 'trifcc2in2'):


		A=triangularfcc2in2 (xini,xfin,yini,yfin,zini,zfin,Npix,Npiy,Ntam,Dx1,Dy1,Dz1,radiusc,radiuscore,radiusshell,A,corenumber,shellnumber,cnum,dl,x,y,z,shape,base,width,height,base2,width2,cs,radiusx1,radiusy1,radiusz1,radiusx2,radiusy2,radiusz2,radiusxcs,radiusycs,radiuszcs)




#	pygeneral (m, ...)
#	!!!

#	print(A)




	return A,Dx1,Dy1,Dz1,zini,zfin,yini,yfin,xini,xfin

###########################################################################################################################


###########################################################################################################################














