/*
 * Declaration of the domain class to manage the physical domain,
 * read the input files and store all the input parameters.
 * Jonas D. De Basabe (jonas@cicese.mx)
 * File:   cdomain.{cpp,h}
 */

#include "macros.h"
#include "cdomain.h"
#include "calgebra.h"
#include "celemesh.h"
#include <stdio.h>

#ifdef ENABLE_EXODUS
#include <exodusII.h>
#endif

extern char LOGFILE[];

char PARFILE[]= INPUTDIR "default.inp";

// Constructor
cDomain::cDomain(char *parfile)//: EPSILON( 10.0*std::numeric_limits<float>::epsilon() )
{// Initialize the properties with the input files
    Init();
    if( parfile==NULL )
        parfile= PARFILE;

    if( ReadParFile(parfile) )
    {
    	//CheckBlocks();
        PhysRec();// Set the physical domain
        //printf("MeshFac=%i\n",MeshFac);
        if( MeshFac>1 ) // Multiply the number of elements and nodes by MeshFac in each direction
        {
            nElem*=MeshFac;
            n*= MeshFac;
            for( int b=0; b<nBlocks; ++b )
                blk[b].nElem*= MeshFac;
        }
    }
    else
    {// There was an error opening the parameters file
        MPIPRINTF("\n\tcDomain ERROR: Unable to read parameters File <%s>\n\n", parfile);
        exit(-1);
    }
}

// Destructor
cDomain::~cDomain(void)
{
	for(int i=0; i<nBlocks; ++i)
		delete blk[i].C;
    delete []vpfile;
    delete []vsfile;
    delete []rhofile;
    delete []meshfile;
    delete []logfile;
    delete []outdir;
}

// Initialize Properties
void cDomain::Init(void)
{
    ExeMode= EXMOD_SET;  // do basic initialization by default
    DispInfo= 1; // Display basic information by default
    vpfile= vsfile= rhofile= meshfile= logfile= outdir= NULL;
    bElast= bFreeSurf= false;
    bHomog= true;
    bSaveMesh= true;
    RandSeed= 0;
    tRandFrac= RF_NO;
    bFracZone= false;
    bSeismoDXU= false;
    nBlocks= nFracs= 0;
    nTaper= 1;
    n.Init(1);
    nElem.Init(1);
    Method= METHOD_SEM;
    TSMethod= TSM_FD;
    MeshType= MESH_REC;
    MeshFac= 1;
    Boundary= BC_NEU;
    tapalpha= tapbeta= 0.0;

    ElemOrd= 2;
    nDim= 2;
    BasisType= NODALGLL;
    vphom= rhohom= 1.0;
    //Vmax= 0.0;
    vshom= 0.0;
    //del.x= del.y= del.z= 1.0;
    delfrac.x= delfrac.y= delfrac.z= 0.0;
    dt= cfl= T= 0.0;
    nt= 0;

    gd.Type= 0;
    gd.ContTime= true;
    gd.Plot= false;

    // DG
    DG.Type= DGM_SIPG;
    DG.bExactInt= EXACT_INT;
    DG.Penalty= SIGMA;
    DG.bPenVel= false;

    // EG
    EG.Ord= 0;
    EG.Type= NODALGLL;

    // Default Source parameters
    src.nShots= 1;
    src.tFunc= 2;
    src.bSave= false;
    src.PkFq= 15.0;
    src.tWidth= 1.0;
    src.sWidth.x = src.sWidth.y = src.sWidth.z= 1.0;
    src.loc.x= src.loc.y= src.loc.z= 0.0;
    src.nor.x= 1.0;
    src.nor.y= 0.0;
    src.nor.z= 0.0;
    src.Amp= 1.0;
    src.Type= SRC_PTSRC;

    // Defalt Receiver parameters
    rec.nTraces= rec.nSeismo= 0;
    rec.r0.x= rec.r0.y= rec.r0.z= 0.0;
    rec.d.x= rec.d.y= rec.d.z= 1.0;
    Snapsh= SS_NONE;

    // Domain
    PT0.x= PT0.y= PT0.z= 0.0;
    PT1= PT0;
    for(int i=0; i<8; ++i)
        PT[i]= PT0;

    // MPI
    MPI_id = MPI_np= 0;
}

// Display the parameters
void cDomain::Display(void)
{
    try
    {
        const char *TSName[3]= {"FD-2", "RK-4", "TS-4"};
        dVec3D *buff_pt0= new dVec3D[MPI_np],
            *buff_pt1= new dVec3D[MPI_np];
        FILE *hFile= fopen(logfile,"at");
#ifdef ENABLEMPI
        MPI::COMM_WORLD.Gather(&DD.PT0, 3, MPI::DOUBLE, buff_pt0, 3, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Gather(&DD.PT1, 3, MPI::DOUBLE, buff_pt1, 3, MPI::DOUBLE, 0);
#endif
        if( MPI_id==0 )
        {
            if( DispInfo )
            {
                if( outdir!=NULL )
                    printf("\tOutput Folder: %s\n",outdir);
                if( logfile!=NULL )
                    printf("\tLog file: %s\n",logfile);
#ifdef ENABLEMPI
                printf("\tUsing MPI\n");
#endif
#ifdef ENABLE_EXODUS
                printf("\tUsing Exodus II library\n");
#endif
#ifdef ENABLENETCDF
                printf("\tUsing NetCDF library\n");
#endif
#ifdef ENABLEMETIS
                printf("\tUsing Metis library\n");
#endif
                if( nDim==3 )
                    printf("\tPhysical Domain      = Rec( (%g, %g, %g), (%g, %g, %g) )\n",
                        PT[0].x,PT[0].y,PT[0].z,PT[7].x,PT[7].y,PT[7].z);
                else
                    printf("\tPhysical Domain      = Rec( (%g, %g), (%g, %g) )\n",
                        PT[0].x,PT[0].z,PT[3].x,PT[3].z);
    #ifdef ENABLEMPI
                for( int i=0; i<MPI_np; ++i )
                {
                    if( nDim==3 )
                        printf("\tSubdomain [P-%.2i] = Rec( (%g, %g, %g), (%g, %g, %g) )\n", i,
                            buff_pt0[i].x,buff_pt0[i].y,buff_pt0[i].z,buff_pt1[i].x,buff_pt1[i].y,buff_pt1[i].z);
                    else
                        printf("\tSubdomain [P-%.2i] = Rec( (%g, %g), (%g, %g) )\n", i,
                            buff_pt0[i].x,buff_pt0[i].z,buff_pt1[i].x,buff_pt1[i].z);
                }
    #endif
                // Time-stepping parameters
                printf("\tTime stepping method = %s, ",TSName[TSMethod]);
                if( cfl>0.0 )
                    printf("CFL = %g ", cfl);
                if( T>0.0 )
                    printf("T = %g ", T);
                if( nt>0 )
                    printf("Nt = %i ", nt);
                if( dt>0.0 )
                    printf("dt = %g ", dt);
                printf("\n");
                // Display source parameters
                printf("\tTime Function        = Ricker-%i, Peak Frequency = %g\n", src.tFunc, src.PkFq);
                printf("\tSource Type          = %i\n", src.Type);
                if( nDim==3 )
                {
                    printf("\tSource Width (t,x)   = (%g,%g,%g,%g)\n", src.tWidth,
                        src.sWidth.x, src.sWidth.y, src.sWidth.z);
                    printf("\tSource location      = (%g,%g,%g)\n", src.loc.x, src.loc.y, src.loc.z);
                }
                else
                {
                    printf("\tSource Width (t,x)   = (%g,%g,%g)\n", src.tWidth,
                        src.sWidth.x, src.sWidth.z);
                    printf("\tSource location      = (%g,%g)\n", src.loc.x, src.loc.z);
                }
                // Output parameters
                printf("\tNumber of Traces     = %i\n",rec.nTraces);
                printf("\tNumber of Seismograms= %i\n",rec.nSeismo);
                printf("\tSeismograms type = ");
                if( bSeismoDXU )
                    printf("DX (DAS)\n");
                else
                    printf("U (DISP)\n");
                if( nDim==3 )
                    for( int i=0; i<rec.nSeismo; ++i )
                    {
                        printf("\tS%i (%5g, %5g, %5g) ", i+1, rec.Seismo[i].x, rec.Seismo[i].y, rec.Seismo[i].z);
                        if( (i%2)==1 || i==(rec.nSeismo-1) )
                            printf("\n");
                        else
                            printf("\t");
                    }
                else
                    for( int i=0; i<rec.nSeismo; ++i )
                    {
                        printf("\tS%i (%5g, %5g) ", i+1, rec.Seismo[i].x, rec.Seismo[i].z);
                        if( (i%2)==1 || i==(rec.nSeismo-1) )
                            printf("\n");
                        else
                            printf("\t");
                    }
                printf("\n");
                // Boundary conditions
                printf("\tBoundary Type        = %i\n",Boundary);
                if( Boundary==BC_TAPER )
                    printf("\t\tTaper Width = %i, alpha = %g, beta = %g\n",
                        nTaper,tapalpha,tapbeta);

                // Layers
                printf("\tNumber of Blocks     = %i\n",nBlocks);
                if( nDim==3 )
                    for( int i=0; i<nBlocks; ++i )
                    {
                        printf("\tBlock(%i) [(1)(%g,%g,%g), (4)(%g,%g,%g), (8)(%g,%g,%g)]",i,
                            blk[i].pt[0].x,blk[i].pt[0].y,blk[i].pt[0].z,
                            blk[i].pt[3].x,blk[i].pt[3].y,blk[i].pt[3].z,
                            blk[i].pt[7].x,blk[i].pt[7].y,blk[i].pt[7].z);
                    	if(blk[i].C)
                    		printf(" Anisotropic\n");
                    	else
                    		printf("\n");
                    }
                else
                    for( int i=0; i<nBlocks; ++i )
                    {
                        printf("\tBlock (%i) [ (%g,%g), (%g,%g) ]",i,
                            blk[i].pt[0].x,blk[i].pt[0].z,
                            blk[i].pt[3].x,blk[i].pt[3].z);
                    	if(blk[i].C)
                    		printf(" Anisotropic\n");
                    	else
                    		printf("\n");
                    }
                printf("\n");
            }

            // Write parameters to the LOG file
            if( hFile )
            {
                if( outdir!=NULL )
                    fprintf(hFile, "\tOutput Folder: %s\n",outdir);
                if( logfile!=NULL )
                    fprintf(hFile, "\tLog file: %s\n",logfile);
                if( nDim==3 )
                    fprintf(hFile, "\tPhysical Domain      = Rec( (%g, %g, %g), (%g, %g, %g) )\n",
                        PT0.x,PT0.y,PT0.z,PT1.x,PT1.y,PT1.z);
                else
                    fprintf(hFile, "\tPhysical Domain      = Rec( (%g, %g), (%g, %g) )\n",
                        PT0.x,PT0.z,PT1.x,PT1.z);
    #ifdef ENABLEMPI
                for( int i=0; i<MPI_np; ++i )
                {
                    if( nDim==3 )
                        fprintf(hFile, "\tSubdomain [P-%.2i] = Rec( (%g, %g, %g), (%g, %g, %g) )\n", i,
                            buff_pt0[i].x,buff_pt0[i].y,buff_pt0[i].z,buff_pt1[i].x,buff_pt1[i].y,buff_pt1[i].z);
                    else
                        fprintf(hFile, "\tSubdomain [P-%.2i] = Rec( (%g, %g), (%g, %g) )\n", i,
                            buff_pt0[i].x,buff_pt0[i].z,buff_pt1[i].x,buff_pt1[i].z);
                }
    #endif
                fprintf(hFile,"\tSimulation time      = %g, %i x %g\n",T,nt, dt);
                fprintf(hFile,"\tTime stepping method = %s\n",TSName[TSMethod]);
                // Display source parameters
                fprintf(hFile,"\tTime Function        = Ricker-%i, Peak Frequency = %g\n", src.tFunc, src.PkFq);
                if( nDim==3 )
                {
                    fprintf(hFile,"\tSource Width (t,x)   = (%g,%g,%g,%g)\n", src.tWidth,
                        src.sWidth.x, src.sWidth.y, src.sWidth.z);
                    fprintf(hFile,"\tSource location      = (%g,%g,%g)\n", src.loc.x, src.loc.y, src.loc.z);
                }
                else
                {
                    fprintf(hFile,"\tSource Width (t,x)   = (%g,%g,%g)\n", src.tWidth,
                        src.sWidth.x, src.sWidth.z);
                    fprintf(hFile,"\tSource location      = (%g,%g)\n", src.loc.x, src.loc.z);
                }
                // Output parameters
                fprintf(hFile,"\tNumber of Traces     = %i\n",rec.nTraces);
                fprintf(hFile,"\tNumber of Seismograms= %i\n",rec.nSeismo);
                if( nDim==3 )
                    for( int i=0; i<rec.nSeismo; ++i )
                    {
                        fprintf(hFile,"\tS%i (%5g, %5g, %5g)\n", i+1, rec.Seismo[i].x, rec.Seismo[i].y, rec.Seismo[i].z);
                        if( (i%2)==1 || i==(rec.nSeismo-1) )
                            fprintf(hFile,"\n");
                        else
                            fprintf(hFile,"\t");
                    }
                else
                    for( int i=0; i<rec.nSeismo; ++i )
                    {
                        fprintf(hFile,"\tS%i (%5g, %5g) ", i+1, rec.Seismo[i].x, rec.Seismo[i].z);
                        if( (i%2)==1 || i==(rec.nSeismo-1) )
                            fprintf(hFile,"\n");
                        else
                            fprintf(hFile,"\t");
                    }
                fprintf(hFile,"\tBoundary Type        = %i\n",Boundary);
                if( Boundary==BC_TAPER )
                    fprintf(hFile,"\t\tTaper Width = %i, alpha = %g, beta = %g\n",
                        nTaper,tapalpha,tapbeta);
                fprintf(hFile,"\tNumber of Blocks     = %i\n",nBlocks);
                if( nDim==3 )
                    for( int i=0; i<nBlocks; ++i )
                    {
                        fprintf(hFile,"\tBlock (%i) [ (1)(%g,%g,%g), (4)(%g,%g,%g), (8)(%g,%g,%g) ]\n",i,
                            blk[i].pt[0].x,blk[i].pt[0].y,blk[i].pt[0].z,
                            blk[i].pt[3].x,blk[i].pt[3].y,blk[i].pt[3].z,
                            blk[i].pt[7].x,blk[i].pt[7].y,blk[i].pt[7].z);
                    }
                else
                    for( int i=0; i<nBlocks; ++i )
                    {
                        fprintf(hFile,"\tBlock (%i) [ (%g,%g), (%g,%g) ]\n",i,
                            blk[i].pt[0].x,blk[i].pt[0].z,
                            blk[i].pt[3].x,blk[i].pt[3].z);
                    }
                fclose(hFile);
            }
        }
        delete[]buff_pt0;
        delete[]buff_pt1;
    }
    catch(int err)
    {
        printf("MPI error %i!\n",err);
    }
}

// Get the domain
bool cDomain::PhysRec(void)
{
    int nn= (nDim==2) ? 4 : 8;
    if( nBlocks<1 )
    {
        MPIPRINTF("cDomain::PhysRec error - Inconsistent block definition");
        return false;
    }
    PT0= blk[0].pt[0];
    PT1= blk[0].pt[0];
    for(int i=0; i<nn; ++i)
        PT[i]= blk[0].pt[i];

    for(int b=0; b<nBlocks; ++b)
    {
        for(int i=0; i<nn; ++i)
        {
            if( PT0.x>blk[b].pt[i].x )
                PT0.x= blk[b].pt[i].x;
            if( PT0.y>blk[b].pt[i].y )
                PT0.y= blk[b].pt[i].y;
            if( PT0.z>blk[b].pt[i].z )
                PT0.z= blk[b].pt[i].z;

            if( PT1.x<blk[b].pt[i].x )
                PT1.x= blk[b].pt[i].x;
            if( PT1.y<blk[b].pt[i].y )
                PT1.y= blk[b].pt[i].y;
            if( PT1.z<blk[b].pt[i].z )
                PT1.z= blk[b].pt[i].z;
        }
        if( nDim==2 && blk[b].Type==BT_LAYER )
        {
            if( PT[0].z>blk[b].pt[0].z )
            {
                PT[0]= blk[b].pt[0];
                PT[1]= blk[b].pt[1];
            }
            if( PT[2].z<blk[b].pt[2].z )
            {
                PT[2]= blk[b].pt[2];
                PT[3]= blk[b].pt[3];
            }
        }
        if( nDim==3 && blk[b].Type==BT_LAYER )
        {
            if( PT[0].z>blk[b].pt[0].z )
            {
                for(int i=0; i<4; ++i)
                    PT[i]= blk[b].pt[i];
            }
            if( PT[2].z<blk[b].pt[2].z )
            {
                for(int i=4; i<8; ++i)
                    PT[i]= blk[b].pt[i];
            }
        }
    }
    return true;
}

double cDomain::GetTimeInc(double h, double v)
{
    try
    {   // First check if any of the time-stepping parameters are defined
        if( dt<=0.0 && nt<=0 && T<=0.0 && cfl<=0.0 )
            throw(1);
        // use dt and nt
        if( dt>0.0 && nt>0 )
            T= dt*nt;
        else if( T>0.0 && nt>0 )// use T and nt
        {
            dt= T/double(nt);
            if( dt*nt < T )
                dt= T/double(nt - 1);
            T= dt*nt;
        }
        else if( T>0.0 && dt>0.0 )// use dt and T
        {
            nt= int(T/dt);
            if( dt*nt < T )
                ++nt;
            T= dt*nt;
        }
        else if( cfl>0.0 && nt>0 )// use cfl and nt
        {
            if( h<=0.0 || v<= 0.0 )
                throw(3);
            dt= cfl*h/v;
            T= dt*nt;
        }
        else if( cfl>0.0 && T>0.0 )// use cfl and T
        {
            if( h<=0.0 || v<= 0.0 )
                throw(3);
            dt= cfl*h/v;
            nt= int(T/dt);
            if( dt*nt < T )
                ++nt;
            T= dt*nt;
        }
        else if( cfl>0.0 && dt>0.0 )// use cfl and dt (not enough info)
            throw(2);
        // Check if all is defined
        if( T<=0.0 || nt<=0 || dt<=0.0 )
            throw(2);
    }
    catch(...)
    {
        MPIPRINTF("\n\tcDomain ERROR: Unable to determine the time step, check input file\n\n");
        return -1.0;
    }
    return dt;
}

// Determine the smallest space increment
void cDomain::GetDelta(double &h, int &b)
{
    double h1, h2, h3;
    h= std::numeric_limits<double>::max();
    printf(" HH=%g \n\n",h);
    for( int i=0; i<nBlocks; ++i )
    {
        if( blk[b].Type!=BT_LAYER )
            continue;

        if( nDim==2 )
        {
            h1= (blk[i].pt[3].x - blk[i].pt[0].x)/blk[i].nElem.x;
            h3= (blk[i].pt[3].z - blk[i].pt[0].z)/blk[i].nElem.z;
            if( h3<h1 )
                h1= h3;
        }
        else
        {
            h1= (blk[i].pt[7].x - blk[i].pt[0].x)/blk[i].nElem.x;
            h2= (blk[i].pt[7].y - blk[i].pt[0].y)/blk[i].nElem.y;
            h3= (blk[i].pt[7].z - blk[i].pt[0].z)/blk[i].nElem.z;
            if( h2<h1 && h2<h3 )
                h1= h2;
            else if( h3<h1 && h3<h2 )
                h1= h3;
        }
        if( h1<h )
        {
            h= h1;
            b= i;
        }
    }
}

// Determine the smallest space increment
void cDomain::GetDelta(double &h, double &v)
{
	double Vmax, h1, h2, h3, t1, tt=std::numeric_limits<double>::max();
    h= tt;
    printf(" HH=%g \n\n",h);

    for( int i=0; i<nBlocks; ++i )
    {// Skip blocks not defined as layers
        if( blk[i].Type!=BT_LAYER )
            continue;

        if( nDim==2 )
        {
            h1= (blk[i].pt[3].x - blk[i].pt[0].x)/blk[i].nElem.x;
            h3= (blk[i].pt[3].z - blk[i].pt[0].z)/blk[i].nElem.z;
            if( h3<h1 )
                h1= h3;
        }
        else
        {
            h1= (blk[i].pt[7].x - blk[i].pt[0].x)/blk[i].nElem.x;
            h2= (blk[i].pt[7].y - blk[i].pt[0].y)/blk[i].nElem.y;
            h3= (blk[i].pt[7].z - blk[i].pt[0].z)/blk[i].nElem.z;
            if( h2<h1 && h2<h3 )
                h1= h2;
            else if( h3<h1 && h3<h2 )
                h1= h3;
        }
        Vmax= blk[i].Vmax;
        Vmax= (is_small(Vmax))? blk[i].Vp : blk[i].Vmax;
        t1= h1/Vmax;
        if( t1<tt )
        {
            h= h1;
            tt= t1;
            v= Vmax;
        }
    }
}

// Get the Maximum Velocity
double cDomain::GetMaxV(void)
{
    double vmax=0.0;
    for( int b=0; b<nBlocks; ++b )
    {
        if( blk[b].Rho!=0.0 )
        {
            if( vmax<blk[b].Vmax )
                vmax= blk[b].Vmax;
            if( vmax<blk[b].Vp )
                vmax= blk[b].Vp;
        }
    }
    if( !Vp.IsEmpty() )
    {
        double vv= Vp.Max();
        if( vmax<vv )
            vmax= vv;
    }
    return vmax;
}

// Get the Maximum Velocity
double cDomain::GetMaxV(int b)
{
    double vmax=0.0;
    if( blk[b].Rho!=0.0 )
    {
        if( vmax<blk[b].Vmax )
            vmax= blk[b].Vmax;
        if( vmax<blk[b].Vp )
            vmax= blk[b].Vp;
    }

    return vmax;
}

// Compute the Anisotropic elastic coefficients using Schoenber-Muir theory
bool cDomain::SetFracEqMedia(void)
{
	double lam, mu, l2m, x1, x2, h=delfrac.Norm(),
		Zn= frac[0].Zn/h, Zt= frac[0].Zt/h;
	cElasTensor *C= NULL;

	if( nBlocks<0 || nFracs<0 )
		return false;
	CalcMediaLame(lam,mu,l2m,0);
	if( blk[0].C==NULL )
		blk[0].C= new cElasTensor();
	C= blk[0].C;
    if( C->c(1,1)==0.0 )
    {
        C->SetLame(lam,mu);
        //C->Display();
        //printf("\n");
        l2m= lam+2.0*mu;
        // Increments on CTT
        x1= Zn/(Zn*l2m+1.0);
        //MPIPRINTF("(l2m=%g)(Zn=%g)(x1=%g)", l2m, Zn, x1);
        x2= Zt/(Zt*mu +1.0);
        C->Add(1,1, -x1*sqr(lam) );
        C->Add(1,2, -x1*sqr(lam) );
        C->Add(2,2, -x1*sqr(lam) );
        // Increments on CNN
        x1*= l2m;
        //MPIPRINTF("(x1=%g)",x1);
        x2*= mu;
        C->Add(3,3, -x1*l2m );
        C->Add(4,4, -x2*mu );
        C->Add(5,5, -x2*mu );
        // Increments on CTN
        C->Add(1,3, -lam*x1);
        C->Add(2,3, -lam*x1);
        if( DispInfo>0 && MPI_id==0 )
        {
            MPIPRINTF("\tSetting up the Anisotropic Equivalent Media...\n");
            MPIPRINTF("\th=%g, Zn=%g, Zt= %g, lambda= %g, mu= %g\n\t", h, Zn, Zt, lam, mu);
            C->Display();
        }
    }
	return true;
}
// Returns the number of time steps between saves
int cDomain::GetNTSave(void)
{
    if( dt<= 0.0 )
        return -1;
    if( nTSave<=0 )
        nTSave= int(fTSave/dt);
    return nTSave;
}
// Read the input file
bool cDomain::GetXBlocks(cVector<int> &XBlks)
{
    int i, b, nXBlks=0;

    for( b=0, i=0; b<nBlocks; ++b) // How many blocks there are in the domain
    {
        if (blk[b].Type==BT_BLOCK)
            i++;
    }
    nXBlks= i;

    if (i<=0)
    {// No hanging blocks
        XBlks.Empty();
        return false;
    }

    XBlks.Init(nXBlks);

    for ( b=0, i=0; b<nBlocks; ++b)
    {
        if(blk[b].Type==BT_BLOCK && i<nXBlks) // Save block in the vector TotalBlocks
        {
            XBlks[i]=b;
            ++i;
        }
    }
    if( DispInfo>2 )
    {
        printf("\n\tHanging Blocks: ");
        XBlks.Display();
    }
    return true;
}
// Read the input file
void cDomain::InitBlocks(void)
{
	blk.Init(nBlocks);
	for(int i,b=0; b<nBlocks; ++b)
	{
		blk[b].bHomog= true;
		blk[b].nElem.Init(0);
		blk[b].Vp= blk[b].Vs= blk[b].Rho= blk[b].Vmax= blk[b].Penalty= 0.0;
		blk[b].C= NULL;
        blk[b].h.Init(0.0);
        blk[b].nFrac= 0;
		blk[b].tFrac= 0;
        blk[b].FracZt= blk[b].FracZn= 0.0;
        blk[b].Type= BT_LAYER;
        for(i=0; i<8; ++i)
        	blk[b].pt[i].Init(0.0);
	}
}

// Check the blocks, return false if there are errors
bool cDomain::CheckBlocks(void)
{
    int b;

    if( nBlocks<1 )
        return false;
    nElem.Init(0);
    // Test #1: check if all the parameters where read
    if( nDim==2 )
        for( b=0; b<nBlocks; ++b )
        {
        	//blk[b].C->Display();
        	if( blk[b].h.x>0.0 && blk[b].h.z>0.0 )
        	{
        		blk[b].pt[1].x= blk[b].pt[0].x + blk[b].h.x;
        		blk[b].pt[1].z= blk[b].pt[0].z;
        		blk[b].pt[2].x= blk[b].pt[0].x;
        		blk[b].pt[2].z= blk[b].pt[0].z + blk[b].h.z;
        		blk[b].pt[3].x= blk[b].pt[0].x + blk[b].h.x;
        		blk[b].pt[3].z= blk[b].pt[0].z + blk[b].h.z;
        	}
            if( blk[b].nElem.x<1 || blk[b].nElem.z<1 )
                return false;
            if( blk[b].Type==BT_LAYER )
            {
                if( nElem.x<blk[b].nElem.x )
                    nElem.x= blk[b].nElem.x;
                nElem.z+= blk[b].nElem.z;
            }
        }
    else
        for( b=0; b<nBlocks; ++b )
        {
        	if( blk[b].h.x>0.0 && blk[b].h.y>0.0 && blk[b].h.z>0.0 )
        	{
        		blk[b].pt[1].x= blk[b].pt[0].x + blk[b].h.x;
        		blk[b].pt[1].y= blk[b].pt[0].y;
        		blk[b].pt[1].z= blk[b].pt[0].z;

        		blk[b].pt[2].x= blk[b].pt[0].x;
        		blk[b].pt[2].y= blk[b].pt[0].y + blk[b].h.y;
        		blk[b].pt[2].z= blk[b].pt[0].z;

        		blk[b].pt[3].x= blk[b].pt[0].x + blk[b].h.x;
        		blk[b].pt[3].y= blk[b].pt[0].y + blk[b].h.y;
        		blk[b].pt[3].z= blk[b].pt[0].z;

        		blk[b].pt[4].x= blk[b].pt[0].x;
        		blk[b].pt[4].y= blk[b].pt[0].y;
        		blk[b].pt[4].z= blk[b].pt[0].z + blk[b].h.z;

        		blk[b].pt[5].x= blk[b].pt[0].x + blk[b].h.x;
        		blk[b].pt[5].y= blk[b].pt[0].y;
        		blk[b].pt[5].z= blk[b].pt[0].z + blk[b].h.z;

        		blk[b].pt[6].x= blk[b].pt[0].x;
        		blk[b].pt[6].y= blk[b].pt[0].y + blk[b].h.y;
        		blk[b].pt[6].z= blk[b].pt[0].z + blk[b].h.z;

        		blk[b].pt[7].x= blk[b].pt[0].x + blk[b].h.x;
        		blk[b].pt[7].y= blk[b].pt[0].y + blk[b].h.y;
        		blk[b].pt[7].z= blk[b].pt[0].z + blk[b].h.z;
        	}
            if( blk[b].nElem.x<1 || blk[b].nElem.y<1 || blk[b].nElem.z<1 )
                return false;
            if( blk[b].Type==BT_LAYER )
            {
                if( nElem.x<blk[b].nElem.x )
                    nElem.x= blk[b].nElem.x;
                if( nElem.y<blk[b].nElem.y )
                    nElem.y= blk[b].nElem.y;
                nElem.z+= blk[b].nElem.z;
            }
        }
    // Test #2: Check if there are gaps or overlaps between the blocks
    //for( b=1; b<nBlocks; ++b )
    //{}
    return true;
}

// Read the input file
bool cDomain::ReadParFile(char *parfile)
{
    ifstream fs;
    char buff[MAXBUFF];
    int aux;

    MPIPRINTF("\tInput file: %s\n",parfile);
    try
    {
        fs.open(parfile, ios::in);
        if( !fs.is_open() )
            throw(1);
        bHomog= true;
        while( !fs.eof() )
        {
            fs.getline(buff,MAXBUFF);
            if( buff[0] != '!' &&  buff[0] != '%' &&  buff[0] != '#' &&  buff[0] != ' ' &&  buff[0] != '\0'
                &&  buff[0] != '\n' &&  buff[0] != '\r' &&  buff[0] != '\t' )
                    ExtractParam(buff);
        }
        fs.close();
        if( Method==METHOD_EG || Method==METHOD_HEG )
        {
            BasisType= NODALENRICHED;
            //ElemOrd= 2; // Only linear basis for enriched Galerkin
        }
        if((Model==MODEL_FRAC || Model==MODEL_ANISOFRAC) && nFracs>1
                && (delfrac.x>0.0 || delfrac.y>0.0 || delfrac.z>0.0) )
        {
            //MPIPRINTF("\tParallel fractures = %i, delta = (%f,%f,%f)\n",
            //        nFracs, delfrac.x, delfrac.y, delfrac.z);
            for(int f=1; f<nFracs; ++f)
            {
                frac[f].Zt= frac[0].Zt;
                frac[f].Zn= frac[0].Zn;
                for(int i=0; i<4; ++i)
                {
                    frac[f].pt[i].x= frac[f-1].pt[i].x + delfrac.x;
                    frac[f].pt[i].y= frac[f-1].pt[i].y + delfrac.y;
                    frac[f].pt[i].z= frac[f-1].pt[i].z + delfrac.z;
                }
            }
        }
        if( !bHomog )
        {
            Vp.Init(n.z, n.x, 0.0);
            Vs.Init(n.z, n.x, 0.0);
            Rho.Init(n.z, n.x, 0.0);
            if( vpfile==NULL )
                Vp= vphom;
            else
                ReadMediaPar(vpfile, Vp);
            if( vsfile==NULL )
                Vs= vshom;
            else
                ReadMediaPar(vsfile, Vs);
            if( rhofile==NULL )
                Rho= rhohom;
            else
                ReadMediaPar(rhofile, Rho);
        }
        if( outdir==NULL )
        {
            //outdir= (char*) calloc(strlen(DEFOUTDIR)+1,sizeof(char));
            aux= strlen(DEFOUTDIR);
            outdir= new char[aux+1];
            strncpy(outdir,DEFOUTDIR,aux);
            outdir[aux]= '\0';
        }
        if( strlen(LOGFILE)>1 )
        {
            delete[]logfile;
            //logfile= (char*) calloc(strlen(LOGFILE)+1,sizeof(char));
            aux= strlen(LOGFILE);
            logfile= new char[aux+1];
            strncpy(logfile,LOGFILE,aux);
            logfile[aux]= '\0';
        }
        else if( logfile!=NULL && strlen(logfile)>1 )
        {
                strncpy(LOGFILE,logfile,strlen(logfile)+1);
        }
        if( CheckBlocks()==false && meshfile!=NULL )
        	throw(2);
        //if( meshfile!=NULL )
        //    ReadMesh(meshfile);
        return true;
    }
    catch(int err)
    {
    	if( err==1 )
    	{
    		fs.close();
    		MPIPRINTF("Error opening input file %s\n",parfile);
    	}
    	else if( err==2 )
        {
            MPIPRINTF("\nInconsistent or incomplete block definition, check input file %s\n", parfile);
            return false;
        }
    }
    return false;
}

// Extract the information from the input lines
void cDomain::ExtractParam(char *buff)
{
    char tag[MAXBUFF], *value;
    int b, aux;
    //double daux, dauxv[8];

    sscanf(buff,"%s",tag);
    value= strchr(buff,'=');
    value++;

    if( strcmp(tag,"SOURCE")==0 )
    {
    	sscanf(value,"%s",tag);
        ReadParFile(tag);
    }
    else if( strcmp(tag,"DISPLAY")==0 )
    {
        sscanf(value,"%i", &aux);
        DispInfo= aux;
    }
    else if( strcmp(tag,"EXEMODE")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
            case EXMOD_SET:
                ExeMode= EXMOD_SET;
                break;
            case EXMOD_MESH:
                ExeMode= EXMOD_MESH;
                break;
            case EXMOD_CALC:
                ExeMode= EXMOD_CALC;
                break;
            case EXMOD_GRDISP:
                ExeMode= EXMOD_GRDISP;
                break;
            default:
                ExeMode= EXMOD_SET;
                break;
        }
    }
    else if( strcmp(tag,"NTSAVE")==0 )
        sscanf(value,"%i", &nTSave);
    else if( strcmp(tag,"DTSAVE")==0 )
        sscanf(value,"%lg", &fTSave);
    else if( strcmp(tag,"SNAPSHRES")==0 )
        sscanf(value,"%i", &Res.x);
    else if( strcmp(tag,"SNAPSHFMT")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
          case SS_NONE:
            Snapsh= SS_NONE;
            break;
          case SS_BIN:
            Snapsh= SS_BIN;
            break;
          case SS_SLICE:
            Snapsh= SS_SLICE;
            break;
          case SS_CDF:
            Snapsh= SS_CDF;
            break;
          case SS_VTK:
              Snapsh= SS_VTK;
              break;
          case SS_EXO:
              Snapsh= SS_EXO;
              break;
          default:
        	  Snapsh= SS_NONE;
        	  break;
        }
    }
    else if( strcmp(tag,"MODEL")==0 )
    {
        sscanf(value,"%i", &aux);
        switch (aux)
        {
          case MODEL_ACOUS:
            Model= MODEL_ACOUS;
            break;
          case MODEL_ELAST:
            Model= MODEL_ELAST;
            break;
          case MODEL_FRAC:
            Model= MODEL_FRAC;
            break;
          case MODEL_ANISOTROPIC:
            Model= MODEL_ANISOTROPIC;
            break;
          case MODEL_ANISOFRAC:
            Model= MODEL_ANISOFRAC;
            break;
          case MODEL_ELACU:
            Model= MODEL_ELACU;
            break;
          case MODEL_MICROLAYERS:
            Model= MODEL_MICROLAYERS;
            break;
          default:
            Model= MODEL_ACOUS;
            break;
        }
    }
    else if( strcmp(tag,"DIMENSIONS")==0 )
    {
        sscanf(value,"%i", &nDim);
        if( nDim<2 ) nDim= 2;
        if( nDim>3 ) nDim= 3;
    }
    else if( strcmp(tag,"METHOD")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
            case METHOD_SG:
                Method= METHOD_SG;
                break;
            case METHOD_VSFD:
                Method= METHOD_VSFD;
                break;
            case METHOD_SEM:
                Method= METHOD_SEM;
                break;
            case METHOD_DG:
                Method= METHOD_DG;
                break;
            case METHOD_IGA:
                Method= METHOD_IGA;
                break;
            case METHOD_EG:
                Method= METHOD_EG;
                break;
            case METHOD_HEG:
                Method= METHOD_HEG;
                break;
            default:
                Method= METHOD_SEM;
                break;
        }
    }
    else if( strcmp(tag,"TSMETHOD")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
          case TSM_FD:
            TSMethod= TSM_FD;
            break;
          case TSM_RK:
            TSMethod= TSM_RK;
            break;
          case TSM_LW4:
            TSMethod= TSM_LW4;
            break;
          default:
            TSMethod= TSM_FD;
            break;
        }
    }
    else if( strcmp(tag,"BASISTYPE")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
          case NODALEQUI:
            BasisType= NODALEQUI;
            break;
          case NODALGLL:
            BasisType= NODALGLL;
            break;
          case NODALGAUSS:
            BasisType= NODALGAUSS;
            break;
          case MODALLEGEN:
            BasisType= MODALLEGEN;
            break;
          default:
            BasisType= NODALGLL;
            break;
        }
    }
    else if( strcmp(tag,"FEMORD")==0 )
    {
        sscanf(value,"%i", &ElemOrd);
        if( ElemOrd<2 ) ElemOrd= 2;
        if( ElemOrd>MAX_ORD ) ElemOrd= MAX_ORD;
    }
    else if( strcmp(tag,"GDTIME")==0 )
    {
        sscanf(value,"%i", &aux);
        gd.ContTime= (aux==0)? true : false;
    }
    else if( strcmp(tag,"GDPLOT")==0 )
    {
        sscanf(value,"%i", &aux);
        gd.Plot= (aux==0)? false : true;
    }
    else if( strcmp(tag,"GDTYPE")==0 )
    {
        sscanf(value,"%i", &aux);
        gd.Type= aux;
    }
    else if( strcmp(tag,"DGTYPE")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
          case DGM_NIPG:
            DG.Type= DGM_NIPG;
            break;
          case DGM_SIPG:
            DG.Type= DGM_SIPG;
            break;
          case DGM_IIPG:
            DG.Type= DGM_IIPG;
            break;
          case DGM_OBBG:
            DG.Type= DGM_OBBG;
            break;
          default:
            DG.Type= DGM_SIPG;
            break;
        }
    }
    else if( strcmp(tag,"DGPENALTY")==0 )
    {
        sscanf(value,"%lg", &DG.Penalty);
    }
    else if( strcmp(tag,"DGLUMPPING")==0 )
    {
        sscanf(value,"%i", &aux);
        DG.bExactInt= (aux==0)? true : false;
    }
    else if( strcmp(tag,"DGPENVEL")==0 )
    {
        sscanf(value,"%i", &aux);
        DG.bPenVel= (aux!=0)? true : false;
    }
    else if( strcmp(tag,"ENRICHORD")==0 )
    {
        sscanf(value,"%i", &aux);
        EG.Ord= (aux<0)? 0 : aux;
    }
    else if( strcmp(tag,"ENRICHTYPE")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
          case NODALEQUI:
            EG.Type= NODALEQUI;
            break;
          case NODALGLL:
            EG.Type= NODALGLL;
            break;
          case NODALGAUSS:
            EG.Type= NODALGAUSS;
            break;
          case MODALLEGEN:
            EG.Type= MODALLEGEN;
            break;
          default:
            EG.Type= NODALGLL;
            break;
        }

    }
    else if( strcmp(tag,"MESHTYPE")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {//MESH_REC, MESH_QUAD, MESH_CUB, MESH_HEX, MESH_THEX
            case MESH_REC:
                MeshType= MESH_REC;
                break;
            case MESH_QUAD:
                MeshType= MESH_QUAD;
                break;
            case MESH_CUB:
                MeshType= MESH_CUB;
                break;
            case MESH_HEX:
                MeshType= MESH_HEX;
                break;
            case MESH_THEX:
                MeshType= MESH_THEX;
                break;
            default:
                MeshType= MESH_REC;
                break;
        }
    }
    else if( strcmp(tag,"SAVEMESH")==0 )
    {
        sscanf(value,"%i", &aux);
        bSaveMesh= aux ? true : false;
    }
    else if( strcmp(tag,"MESHFILE")==0 )
    {
        sscanf(value,"%s",tag);
        aux= (int) strlen(tag);
        meshfile= new char[aux+1];
        strncpy(meshfile,tag,aux+1);
    }
    else if( strcmp(tag,"FREESURF")==0 )
    {
        sscanf(value,"%i", &aux);
        bFreeSurf= (aux==1)? false : true;
    }
    else if( strcmp(tag,"SAVESRC")==0 )
    {
        sscanf(value,"%i", &aux);
        src.bSave= (aux==0)? false : true;
    }
    else if( strcmp(tag,"NT")==0 )
        sscanf(value,"%i", &nt);
    else if( strcmp(tag,"TMAX")==0 )
        sscanf(value,"%lg", &T);
    else if( strcmp(tag,"DELTAT")==0 )
        sscanf(value,"%lg", &dt);
    else if( strcmp(tag,"CFL")==0 )
        sscanf(value,"%lg", &cfl);
    else if( strcmp(tag,"BC")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
          case BC_NEU:
            Boundary= BC_NEU;
            break;
          case BC_TAPER:
            Boundary= BC_TAPER;
            break;
          case BC_PERIODIC:
            Boundary= BC_PERIODIC;
            break;
          case BC_PARAXABC:
            Boundary= BC_PARAXABC;
            break;
          default:
            Boundary= BC_NEU;
            break;
        }
    }
    else if( strcmp(tag,"TAPERLEN")==0 )
        sscanf(value,"%i", &nTaper);
    else if( strcmp(tag,"TAPALPHA")==0 )
        sscanf(value,"%lg", &tapalpha);
    else if( strcmp(tag,"TAPBETA")==0 )
        sscanf(value,"%lg", &tapbeta);
    else if( strcmp(tag,"SRCTYPE")==0 )
    {
        sscanf(value,"%i", &aux);
        switch(aux)
        {
          case SRC_PTSRC:
            src.Type= SRC_PTSRC;
            break;
          case SRC_PPLANEWV:
            src.Type= SRC_PPLANEWV;
            break;
          case SRC_SPLANEWV:
            src.Type= SRC_SPLANEWV;
            break;
          default:
            src.Type= SRC_PTSRC;
            break;
        }
    }
    else if( strcmp(tag,"SRCFUNC")==0 )
        sscanf(value,"%i", &src.tFunc);
    else if( strcmp(tag,"TWIDTH")==0 )
    {
        sscanf(value,"%lg", &src.tWidth);
    }
    else if( strcmp(tag,"SWIDTH")==0 )
    {
        if( nDim==2 )
            sscanf(value,"%lg %lg", &src.sWidth.x, &src.sWidth.z);
        else
            sscanf(value,"%lg %lg %lg", &src.sWidth.x, &src.sWidth.y, &src.sWidth.z);
    }
    else if( strcmp(tag,"PKFREQ")==0 )
    {
        sscanf(value,"%lg", &src.PkFq);
    }
    else if( strcmp(tag,"SRCAMP")==0 )
    {
        sscanf(value,"%lg", &src.Amp);
    }
    else if( strcmp(tag,"SHOTS")==0 )
        sscanf(value,"%i", &src.nShots);
    else if( strcmp(tag,"SRCLOC")==0 )
    {
        if( nDim==2 )
        {
            sscanf(value,"%lg %lg", &src.loc.x, &src.loc.z);
        }
        else
        {
            sscanf(value,"%lg %lg %lg", &src.loc.x, &src.loc.y, &src.loc.z);
        }
    }
    else if( strcmp(tag,"SRCVECTOR")==0 )
    {
        if( nDim==2 )
        {
            sscanf(value,"%lg %lg", &src.nor.x, &src.nor.z);
        }
        else
        {
            sscanf(value,"%lg %lg %lg", &src.nor.x, &src.nor.y, &src.nor.z);
        }
    }
    else if( strcmp(tag,"SRCTENSOR")==0 )
    {
        if( nDim==2 )
        {
            sscanf(value,"%lg %lg %lg", &src.tensr[0], &src.tensr[1], &src.tensr[2]);
        }
        else
        {
            sscanf(value,"%lg %lg %lg %lg %lg %lg", &src.tensr[0], &src.tensr[1], &src.tensr[2], &src.tensr[3], &src.tensr[4], &src.tensr[5]);
        }
    }
    else if( strcmp(tag,"REC0")==0 )
    {
        if( nDim==2 )
            sscanf(value,"%lg %lg", &rec.r0.x, &rec.r0.z);
        else
            sscanf(value,"%lg %lg %lg", &rec.r0.x, &rec.r0.y, &rec.r0.z);
    }
    else if( strcmp(tag,"RECD")==0 )
    {
        if( nDim==2 )
            sscanf(value,"%lg %lg", &rec.d.x, &rec.d.z);
        else
            sscanf(value,"%lg %lg %lg", &rec.d.x, &rec.d.y, &rec.d.z);
    }
    else if( strcmp(tag,"NREC")==0 )
        sscanf(value,"%i", &rec.nTraces);
    else if( strcmp(tag,"NSEISMO")==0 )
    {
        sscanf(value,"%i", &rec.nSeismo);
        if( rec.nSeismo>0 )
        {
            rec.Seismo.Init(rec.nSeismo);
            for( aux=0; aux<rec.nSeismo; ++aux )
                rec.Seismo[aux].x= rec.Seismo[aux].y= rec.Seismo[aux].z= 0.0;
        }
    }
    else if( strncmp(tag,"SEISMO",6)==0 )
    {
        aux= ch2int(tag[6]);
        if( aux<=rec.nSeismo )
        {
            if( nDim==2 )
                sscanf(value,"%lg %lg", &rec.Seismo[aux-1].x, &rec.Seismo[aux-1].z);
            else
                sscanf(value,"%lg %lg %lg", &rec.Seismo[aux-1].x,
                        &rec.Seismo[aux-1].y, &rec.Seismo[aux-1].z);
        }
    }
    else if( strncmp(tag,"DXUSEIS",7)==0 )
    {
        sscanf(value,"%i", &aux);
        bSeismoDXU= (aux!=0);
    }
    else if( strcmp(tag,"NFRACTURES")==0 )
    {
        sscanf(value,"%i", &nFracs);
        if( nFracs>0 )
        {
            if( frac.Init(nFracs) )
            {
                for( aux=0; aux<nFracs; ++aux )
                    frac[aux].Zn= frac[aux].Zt= 0.0;
            }
            else
                printf("[cDomain::ExtractParam] Error allocating fractures parameters\n");
        }
    }
    else if( strcmp(tag,"TRANDFRAC")==0 )
    {
        sscanf(value,"%i", &aux);
        //tRandFrac= aux<=0 ? RF_NO : aux;
        switch(aux)
        {
            case 1:
                tRandFrac= RF_B;
                break;
            case 2:
                tRandFrac= RF_H;
                break;
            case 3:
                tRandFrac= RF_V;
                break;
            default:
                tRandFrac= RF_NO;
                break;
        }
    }
    else if( strcmp(tag,"FRACZONE")==0 )
    {
        sscanf(value,"%i", &aux);
        bFracZone= (aux==0)? false : true;
    }
    else if( strcmp(tag,"RANDSEED")==0 )
    {
        sscanf(value,"%i", &RandSeed);
    }
    else if( strcmp(tag,"DELTAFRAC")==0 )
    {
        aux= sscanf(value,"%lg %lg %lg", &delfrac.x, &delfrac.y, &delfrac.z);
        if (aux<3)
        {
            delfrac.z= delfrac.y;
            delfrac.y= 0.0;
        }
    }
    else if( strncmp(tag,"ZT_",3)==0 )
    {
        str2int(tag+3,aux);
        if( aux<=nFracs )
        {
            sscanf(value,"%lg", &frac[aux-1].Zt);
        }
    }
    else if( strncmp(tag,"ZN_",3)==0 )
    {
        str2int(tag+3,aux);
        if( aux<=nFracs )
        {
            sscanf(value,"%lg", &frac[aux-1].Zn);
        }
    }
    else if( strncmp(tag,"FRACX_",6)==0 )
    {
        str2int(tag+6,aux);
        if( aux<=nFracs )
        {
            sscanf(value,"%lg %lg %lg %lg",
            		&frac[aux-1].pt[0].x, &frac[aux-1].pt[1].x, &frac[aux-1].pt[2].x, &frac[aux-1].pt[3].x);
        }
    }
    else if( strncmp(tag,"FRACY_",6)==0 )
    {
        str2int(tag+6,aux);
        if( aux<=nFracs )
        {
            sscanf(value,"%lg %lg %lg %lg",
            		&frac[aux-1].pt[0].y, &frac[aux-1].pt[1].y, &frac[aux-1].pt[2].y, &frac[aux-1].pt[3].y);
        }
    }
    else if( strncmp(tag,"FRACZ_",6)==0 )
    {
        str2int(tag+6,aux);
        if( aux<=nFracs )
        {
            sscanf(value,"%lg %lg %lg %lg",
            		&frac[aux-1].pt[0].z, &frac[aux-1].pt[1].z, &frac[aux-1].pt[2].z, &frac[aux-1].pt[3].z);
        }
    }
    else if( strcmp(tag,"FRACZONEX")==0 )
    {
        if( nDim==2 )
            sscanf(value,"%lg %lg %lg %lg",
                &FZ[0].x, &FZ[1].x, &FZ[2].x, &FZ[3].x);
        else // 3D
            sscanf(value,"%lg %lg %lg %lg %lg %lg %lg %lg",
                &FZ[0].x, &FZ[1].x, &FZ[2].x, &FZ[3].x,
                &FZ[4].x, &FZ[5].x, &FZ[6].x, &FZ[7].x);

    }
    else if( strcmp(tag,"FRACZONEY")==0 )
    {
        if( nDim==2 )
            sscanf(value,"%lg %lg %lg %lg",
                &FZ[0].y, &FZ[1].y, &FZ[2].y, &FZ[3].y);
        else // 3D
            sscanf(value,"%lg %lg %lg %lg %lg %lg %lg %lg",
                &FZ[0].y, &FZ[1].y, &FZ[2].y, &FZ[3].y,
                &FZ[4].y, &FZ[5].y, &FZ[6].y, &FZ[7].y);

    }
    else if( strcmp(tag,"FRACZONEZ")==0 )
    {
        if( nDim==2 )
            sscanf(value,"%lg %lg %lg %lg",
                &FZ[0].z, &FZ[1].z, &FZ[2].z, &FZ[3].z);
        else // 3D
            sscanf(value,"%lg %lg %lg %lg %lg %lg %lg %lg",
                &FZ[0].z, &FZ[1].z, &FZ[2].z, &FZ[3].z,
                &FZ[4].z, &FZ[5].z, &FZ[6].z, &FZ[7].z);

    }
    else if( strcmp(tag,"LOGFILE")==0 )
    {
        sscanf(value,"%s", tag);
        aux= (int) strlen(tag);
        logfile= new char[aux+1];
        strncpy(logfile,tag,aux+1);
    }
    else if( strcmp(tag,"OUTDIR")==0 )
    {
        sscanf(value,"%s", tag);
        aux= (int) strlen(tag);
        outdir= new char[aux+1];
        strncpy(outdir,tag,aux+1);
    }
    else if( strcmp(tag,"VPFILE")==0 )
    {
        bHomog= false;
        sscanf(value,"%s", tag);
        aux= (int) strlen(tag);
        vpfile= new char[aux+1];
        strncpy(vpfile,tag,aux+1);
    }
    else if( strcmp(tag,"VSFILE")==0 )
    {
        bHomog= false;
        sscanf(value,"%s",tag);
        aux= (int) strlen(tag);
        vsfile= new char[aux+1];
        strncpy(vsfile,tag,aux+1);
    }
    else if( strcmp(tag,"RHOFILE")==0 )
    {
        bHomog= false;
        sscanf(value,"%s",tag);
        aux= (int) strlen(tag);
        rhofile= new char[aux+1];
        strncpy(rhofile,tag,aux+1);
    }
    else if( strcmp(tag,"NBLOCKS")==0 )
    {
        sscanf(value,"%i", &nBlocks);
        if( nBlocks<1 )
            nBlocks= 1;
        InitBlocks();
    }
    else if( strcmp(tag,"MESHFAC")==0 )
    {
        sscanf(value,"%i", &MeshFac);
        if( MeshFac<1 )
            MeshFac= 1;
    }
    else if( strncmp(tag,"NELEM_",6)==0 )
    {
        aux= ch2int(tag[6]);
        if( nBlocks>0 && aux<=nBlocks )
        {
            if( nDim==2 )
                sscanf(value,"%i %i", &blk[aux-1].nElem.x, &blk[aux-1].nElem.z);
            else
                sscanf(value,"%i %i %i", &blk[aux-1].nElem.x,
                        &blk[aux-1].nElem.y, &blk[aux-1].nElem.z);
        }
    }
    else if( strncmp(tag,"NP",2)==0 )
    {
        if( nDim==2 )
            sscanf(value,"%i %i", &n.x, &n.z);
    }
    else if( strncmp(tag,"X_",2)==0 )
    {
        aux= ch2int(tag[2]);
        if( nBlocks>0 && aux<=nBlocks )
        {
            if( nDim==2 )
			{
				sscanf(value,"%lg %lg %lg %lg",
						&blk[aux-1].pt[0].x, &blk[aux-1].pt[1].x, &blk[aux-1].pt[2].x, &blk[aux-1].pt[3].x);
			}
            else
			{
				sscanf(value,"%lg %lg %lg %lg %lg %lg %lg %lg",
						&blk[aux-1].pt[0].x, &blk[aux-1].pt[1].x, &blk[aux-1].pt[2].x, &blk[aux-1].pt[3].x,
						&blk[aux-1].pt[4].x, &blk[aux-1].pt[5].x, &blk[aux-1].pt[6].x, &blk[aux-1].pt[7].x);
			}
        }
    }
    else if( strncmp(tag,"Y_",2)==0 )
    {
        aux= ch2int(tag[2]);
        if( nBlocks>0 && aux<=nBlocks )
        {
            if( nDim==2 )
			{
				sscanf(value,"%lg %lg %lg %lg",
						&blk[aux-1].pt[0].y, &blk[aux-1].pt[1].y, &blk[aux-1].pt[2].y, &blk[aux-1].pt[3].y);
			}
            else
			{
				sscanf(value,"%lg %lg %lg %lg %lg %lg %lg %lg",
						&blk[aux-1].pt[0].y, &blk[aux-1].pt[1].y, &blk[aux-1].pt[2].y, &blk[aux-1].pt[3].y,
						&blk[aux-1].pt[4].y, &blk[aux-1].pt[5].y, &blk[aux-1].pt[6].y, &blk[aux-1].pt[7].y);
			}
        }
    }
    else if( strncmp(tag,"Z_",2)==0 )
    {
        aux= ch2int(tag[2]);
        if( nBlocks>0 && aux<=nBlocks )
        {
            if( nDim==2 )
			{
				sscanf(value,"%lg %lg %lg %lg",
						&blk[aux-1].pt[0].z, &blk[aux-1].pt[1].z, &blk[aux-1].pt[2].z, &blk[aux-1].pt[3].z);
			}
            else
			{
				sscanf(value,"%lg %lg %lg %lg %lg %lg %lg %lg",
						&blk[aux-1].pt[0].z, &blk[aux-1].pt[1].z, &blk[aux-1].pt[2].z, &blk[aux-1].pt[3].z,
						&blk[aux-1].pt[4].z, &blk[aux-1].pt[5].z, &blk[aux-1].pt[6].z, &blk[aux-1].pt[7].z);
			}
        }
    }
    else if( strncmp(tag,"P_",2)==0 )
    {
        aux= ch2int(tag[2]);
        if( nBlocks>0 && aux<=nBlocks )
        {
            if( nDim==2 )
            {
                sscanf(value,"%lg %lg", &blk[aux-1].pt[0].x, &blk[aux-1].pt[0].z);
            }
            else
            {
                sscanf(value,"%lg %lg %lg", &blk[aux-1].pt[0].x, &blk[aux-1].pt[0].y, &blk[aux-1].pt[0].z);
            }
        }
    }
    else if( strncmp(tag,"H_",2)==0 )
    {
        aux= ch2int(tag[2]);
        if( nBlocks>0 && aux<=nBlocks )
        {
            if( nDim==2 )
            {
                sscanf(value,"%lg %lg", &blk[aux-1].h.x, &blk[aux-1].h.z);
            }
            else
            {
                sscanf(value,"%lg %lg %lg", &blk[aux-1].h.x, &blk[aux-1].h.y, &blk[aux-1].h.z);
            }
        }
    }
    else if( strncmp(tag,"NFRAC_",6)==0 )
    {
        aux= ch2int(tag[6]);
        if( nBlocks>0 && aux<=nBlocks )
            sscanf(value,"%i", &blk[aux-1].nFrac);
    }
    else if( strncmp(tag,"TRANDFRAC_",10)==0 )
    {
        aux= ch2int(tag[10]);
        if( nBlocks>0 && aux<=nBlocks )
            sscanf(value,"%i", &blk[aux-1].tFrac);
    }
    else if( strncmp(tag,"FRACZT_",7)==0 )
    {
        aux= ch2int(tag[7]);
        if( nBlocks>0 && aux<=nBlocks )
            sscanf(value,"%lg", &blk[aux-1].FracZt);
    }
    else if( strncmp(tag,"FRACZN_",7)==0 )
    {
        aux= ch2int(tag[7]);
        if( nBlocks>0 && aux<=nBlocks )
            sscanf(value,"%lg", &blk[aux-1].FracZn);
    }
    else if( strncmp(tag,"VMAX_",5)==0 )
    {
        aux= ch2int(tag[5]);
        if( nBlocks>0 && aux<=nBlocks )
            sscanf(value,"%lg", &blk[aux-1].Vmax);
    }
    else if( strncmp(tag,"VP_",3)==0 )
    {
        sscanf(value,"%lg", &vphom);
        aux= ch2int(tag[3]);
        if( nBlocks>0 && aux<=nBlocks )
            blk[aux-1].Vp= vphom;
    }
    else if( strncmp(tag,"VS_",3)==0 )
    {
        sscanf(value,"%lg", &vshom);
        aux= ch2int(tag[3]);
        if( nBlocks>0 && aux<=nBlocks )
            blk[aux-1].Vs= vshom;
    }
    else if( strncmp(tag,"RHO_",4)==0 )
    {
        sscanf(value,"%lg", &rhohom);
        aux= ch2int(tag[4]);
        if( nBlocks>0 && aux<=nBlocks )
            blk[aux-1].Rho= rhohom;
    }
    else if( strncmp(tag,"CIJ_",4)==0 )
    {
        aux= ch2int(tag[4]);
        if( nBlocks>0 && aux<=nBlocks )
        {
        	double daux[13];
        	--aux;
        	blk[aux].C= new cElasTensor();
        	if( blk[aux].C )
        	{
				int nn= sscanf(value,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
						&daux[0], &daux[1], &daux[2], &daux[3], &daux[4],
						&daux[5], &daux[6], &daux[7], &daux[8],
						&daux[9], &daux[10], &daux[11], &daux[12]);
				if( nn<13 )
					for(int i=nn; i<13; ++i)
						daux[i]= 0.0;
				blk[aux].C->Init(daux);
				//printf("blk=%i ", aux);
				//blk[aux].C->Display();
        	}
        	else
        		MPIPRINTF("Error reading anisotropic elastic parameters\n");
        }
    }
    else if( strncmp(tag,"BLOCKTYPE_",10)==0 )
    {
        b= ch2int(tag[10]);
        if( nBlocks>0 && b<=nBlocks )
        {
            sscanf(value,"%i", &aux);
            switch(aux)
            {
                case BT_LAYER:
                    blk[b-1].Type= BT_LAYER;
                    break;
                case BT_BLOCK:
                default:
                    blk[b-1].Type= BT_BLOCK;
                    break;
            }
        }
    }
    else if( strncmp(tag,"PENALTY_",8)==0 )
    {
        aux= ch2int(tag[8]);
        if( nBlocks>0 && aux<=nBlocks )
            sscanf(value,"%lg", &blk[aux-1].Penalty);
    }
    else MPIPRINTF("\nUnknown tag <%s> in Input file\n",tag);
}

// Read the Media Parameters
bool cDomain::ReadMediaPar(char *szfile, cMat2<float> &a)
{
    char *ext= strrchr(szfile,'.');

    if( strcmp(ext,".txt")==0 )
        return ReadMediaParTxt(szfile, a);
    else
        return ReadMediaParBin(szfile, a);
}

// Read the Media Parameters
bool cDomain::ReadMediaParBin(char *szfile, cMat2<float> &a)
{
    const int FLTSZ= sizeof(float);
    ifstream fs;
    int i,j;
    float aux;
    char fname[MAXBUFF];

    try
    {
        sprintf(fname,INPUTDIR"%s",szfile);
        fs.open(fname, ios::in|ios::binary);
        if( !fs.is_open() )
            throw 0;

        for( j=0; j<n.x; ++j )
            for( i=0; i<n.z; ++i )
            {
                fs.read((char*) &aux,FLTSZ);
                a[i][j]= aux;
            }
        fs.close();
    }
    catch(...)
    {
        fs.close();
        return false;
    }
    return true;
}

// Read the Media Parameters
bool cDomain::ReadMediaParTxt(char *szfile, cMat2<float> &a)
{
    FILE *hFile;
    float faux;
    int i,j, aux;
    char fname[MAXBUFF];

    sprintf(fname,INPUTDIR"%s",szfile);
    if( (hFile= fopen(fname, "rt"))==0 )
        return false;

    for( j=0; j<n.x; ++j )
        for( i=0; i<n.z; ++i )
        {
            aux= fscanf(hFile,"%f",&faux);
            a[i][j]= (aux>0)? faux : 0.0;
        }
    fclose(hFile);
    return true;
}

// Read the Mesh file
bool cDomain::ReadMesh(cEleMesh &Ele, cMat2<int> &ElemBlk)
{
    char *ext= strrchr(meshfile,'.');

    Ele.Type= MESH_CUB;
    if( strcmp(ext,".txt")==0 )
        return ReadMeshTxt(Ele.Anch,Ele,ElemBlk);
    else
        return ReadMeshBin(Ele.Anch,Ele,ElemBlk);
}

// Read the Mesh file
bool cDomain::ReadMeshTxt(dmat &Nod, cEleMesh &Ele, cMat2<int> &ElemBlk)
{
    FILE *hFile;

    try
    {
        double daux1, daux2;
        int i,j,nNodes,nelem,iaux0,iaux1,iaux2,iaux3;
        char fname[MAXBUFF];

        sprintf(fname,INPUTDIR"%s",meshfile);
        if( (hFile= fopen(fname, "rt"))==0 )
            throw(1);

        j= fscanf(hFile,"%i",&nNodes);
        Nod.Init(2,nNodes,0.0);
        for( i=0; i<nNodes; ++i)
        {
            j= fscanf(hFile,"%lg%lg",&daux1,&daux2);
            if( j<1 ) throw(2);
            Nod[0][i]= daux1;
            Nod[1][i]= daux2;
        }

        j= fscanf(hFile,"%i",&nelem);
        Ele.Init(nelem);
        for( i=0; i<nelem; ++i)
        {
            j= fscanf(hFile,"%i%i%i%i",&iaux0,&iaux1,&iaux2,&iaux3);
            if( j<1 ) throw(2);
            Ele[i][0]= iaux0;
            Ele[i][1]= iaux1;
            Ele[i][2]= iaux2;
            Ele[i][3]= iaux3;
        }
        fclose(hFile);
        ElemBlk.Init(1,2,0);
        ElemBlk[0][0]= 0;
        ElemBlk[0][1]= nelem;
    }
    catch(int ex)
    {
        if( ex==2 )
        {    MPIPRINTF("\nSyntax error in mesh file <%s>\n",meshfile);}
        else
        {    MPIPRINTF("\nError opening mesh file <%s>\n",meshfile);}
        fclose(hFile);
        return false;
    }
    return true;
}

//---------------------------------------------------------------------------
// Read the Mesh file in Exodus II format
#define MINMAXH(p0,p1) {h=p0.Dist(p1); \
    if(h<hmin) hmin= h; \
    if(h>hmax) hmax= h;}
bool cDomain::ReadMeshBin(dmat &Nod, cEleMesh &Ele, cMat2<int> &ElemBlk)
{
#ifdef ENABLE_EXODUS
    int CPU_word_size,IO_word_size, exoid, exoerr,
        num_dim, num_nodes, num_elem, num_elem_blk,
        num_node_sets, num_side_sets, *ElemBlkID,
        nElemBlk,nNodPerElem,nAttr, b, nBlocks;
    char title[MAX_LINE_LENGTH+8], ElemType[MAX_STR_LENGTH+8];
    float version, *Nodes_x, *Nodes_y, *Nodes_z, hmin=1.0e10, hmax=0.0, h;
    fVec3D p[4], q[4];

    try
    {
        // Open the EXODUS II file and read header information
        CPU_word_size = sizeof(float); /* float or double */
        IO_word_size = 0; /* use what is stored in file */
        exoid = ex_open (meshfile, EX_READ, &CPU_word_size, &IO_word_size, &version);
        if( exoid<0 ) throw(1);

        exoerr= ex_get_init (exoid, title, &num_dim, &num_nodes,
                            &num_elem, &num_elem_blk, &num_node_sets, &num_side_sets);
        if( exoerr ) throw(2);
        if( DispInfo )
        {
            MPIPRINTF("\tExodusII File: %s\n\tTitle = %s\n\tDim = %i, Nodes = %i, Elements = %i\n",
                    meshfile, title, num_dim, num_nodes, num_elem);
            MPIPRINTF("\tElement Blocks = %i, Node Sets = %i, Side Sets = %i, Wrd Size = %i\n",
                    num_elem_blk, num_node_sets, num_side_sets, IO_word_size);
        }
        // Read the nodes
        Nod.Init(num_dim,num_nodes,0.0);
        Nodes_x= new float[num_nodes];
        if( Nodes_x==NULL )
            throw(20);
        Nodes_y= new float[num_nodes];
        if( Nodes_y==NULL )
            throw(20);
        if( num_dim>2 )
        {
            Nodes_z= new float[num_nodes];
            if( Nodes_z==NULL )
                throw(20);
        }
        else
            Nodes_z= NULL;
        exoerr= ex_get_coord (exoid,Nodes_x,Nodes_y,Nodes_z);
        if(exoerr) throw(2);
        p[0].x= q[0].x= Nodes_x[0];
        p[0].y= q[0].y= Nodes_y[0];
        p[0].z= q[0].z= (num_dim>2)? Nodes_z[0] : 0.0;
        for( int i=0; i<num_nodes; ++i )
        {
            Nod[0][i]= Nodes_x[i];
            Nod[1][i]= Nodes_y[i];
            // Get the minimum and maximum x
            if( p[0].x>Nodes_x[i] )
                p[0].x= Nodes_x[i];
            if( q[0].x<Nodes_x[i] )
                q[0].x= Nodes_x[i];
            // Get the minimum and maximum y
            if( p[0].y>Nodes_y[i] )
                p[0].y= Nodes_y[i];
            if( q[0].y<Nodes_y[i] )
                q[0].y= Nodes_y[i];
            if( num_dim>2 )
            {
                Nod[2][i]= Nodes_z[i];
                // Get the minimum and maximum z
                if( p[0].z>Nodes_z[i] )
                    p[0].z= Nodes_z[i];
                if( q[0].z<Nodes_z[i] )
                    q[0].z= Nodes_z[i];
            }
        }

        PT0= p[0];
        PT1= q[0];
        DD.PT0= DD.VPT0= PT0;
        DD.PT1= DD.VPT1= PT1;
        if( nDim==2 )
        {
            PT[0]= PT[1]= PT0;
            PT[2]= PT[3]= PT1;
            PT[1].x= PT1.x;
            PT[2].x= PT0.x;
            for( int i=0; i<4; ++i )
                DD.PT[i]= PT[i];
        }
        else // nDim==3
        {
            PT[0]= PT[1]= PT[2]= PT[3]= PT0;
            PT[4]= PT[5]= PT[6]= PT[7]= PT1;
            PT[1].x= PT[3].x= PT1.x;
            PT[2].y= PT[3].y= PT1.y;
            PT[4].x= PT[6].x= PT0.x;
            PT[4].y= PT[5].y= PT0.y;
            for( int i=0; i<8; ++i )
                DD.PT[i]= PT[i];
        }
        // Read the elements
        ElemBlkID= new int[num_elem_blk];
        if( ElemBlkID==NULL )
            throw(20);
        exoerr= ex_get_elem_blk_ids(exoid, ElemBlkID);
        if(exoerr) throw(2);

        // Determine the number of elements in the compatible blocks
        num_elem= 0;
        nBlocks= 0;
        for( int i=0; i<num_elem_blk; ++i )
        {
            exoerr= ex_get_elem_block (exoid, ElemBlkID[i], ElemType, &nElemBlk,
                                       &nNodPerElem, &nAttr);
            if(exoerr) throw(2);
            if( ( num_dim==2 && nNodPerElem==4 ) || ( num_dim==3 && nNodPerElem==8 ) )
            {
            	++nBlocks;
            	num_elem+= nElemBlk;
            }
        }
        if( nBlocks<=0 )
        	throw(4);
        ElemBlk.Init(nBlocks,2,0);
        if( num_dim==2 )
        {
        	nNodPerElem= 4;
        	Ele.InitQuad(num_elem);
        }
        else
        {
        	nNodPerElem= 8;
        	Ele.InitHexa(num_elem);
        }
        b= 0;
        for( int i=0,j=0; i<num_elem_blk; ++i )
        {
            int *ElemConn;
            exoerr= ex_get_elem_block (exoid, ElemBlkID[i], ElemType, &nElemBlk,
                                       &nNodPerElem, &nAttr);
            if(exoerr) throw(2);
            if( DispInfo )
                MPIPRINTF("\tBlock[%2i] : Type = %s, Elements = %i, Nodes/Elem = %i, Attr. = %i\n",
                      i, ElemType,nElemBlk,nNodPerElem,nAttr);
            if( ( num_dim==2 && nNodPerElem!=4 ) || ( num_dim==3 && nNodPerElem!=8 ) )
            {
                MPIPRINTF("\t\t>skipping block, incompatible Element Type\n");
            }
            else
            {
				ElemConn= new int[nElemBlk*nNodPerElem];
                if( ElemConn==NULL )
                    throw(20);
				exoerr= ex_get_elem_conn (exoid, ElemBlkID[i], ElemConn);
				if(exoerr) throw(2);
				ElemBlk[b][0]= j;
				for( int e=0,l=0; e<nElemBlk; ++e,++j)
				{
					Ele[j].blk= b; // assign block number
					//printf("[ %i - %i - %i ] ",e,l,j);
					//for( k=0; k<nNodPerElem; ++k, ++l)
					//    elmnds[k]=  ElemConn[l] - 1;
						//Ele[j][k]= ElemConn[l] - 1;
					if( nNodPerElem==4 )
					{
						Ele[j][3]= ElemConn[l++] - 1;
						Ele[j][2]= ElemConn[l++] - 1;
						Ele[j][0]= ElemConn[l++] - 1;
						Ele[j][1]= ElemConn[l++] - 1;
						// Check the distance between the nodes
						for( int i=0; i<4; ++i )
						{
							p[i].x= Nod[0][Ele[j][i]];
							p[i].y= Nod[1][Ele[j][i]];
							p[i].z= (num_dim>2)? Nod[2][Ele[j][i]] : 0.0;
						}
						// Get the minimum and maximum h
						MINMAXH(p[0],p[1]);
						MINMAXH(p[0],p[2]);
						MINMAXH(p[3],p[1]);
						MINMAXH(p[3],p[2]);
					}
					else if( nNodPerElem==8 )
					{
						Ele[j][0]= ElemConn[l++] - 1;
						Ele[j][1]= ElemConn[l++] - 1;
						Ele[j][3]= ElemConn[l++] - 1;
						Ele[j][2]= ElemConn[l++] - 1;
						Ele[j][4]= ElemConn[l++] - 1;
						Ele[j][5]= ElemConn[l++] - 1;
						Ele[j][7]= ElemConn[l++] - 1;
						Ele[j][6]= ElemConn[l++] - 1;
						// Check the distance between the nodes
						for( int i=0; i<4; ++i )
						{
							p[i].x= Nod[0][Ele[j][i]];
							p[i].y= Nod[1][Ele[j][i]];
							p[i].z= Nod[2][Ele[j][i]];
							q[i].x= Nod[0][Ele[j][i+4]];
							q[i].y= Nod[1][Ele[j][i+4]];
							q[i].z= Nod[2][Ele[j][i+4]];
						}
						// Get the minimum and maximum h
						MINMAXH(p[0],p[1]);
						MINMAXH(p[0],p[2]);
						MINMAXH(p[3],p[1]);
						MINMAXH(p[3],p[2]);
						MINMAXH(q[0],q[1]);
						MINMAXH(q[0],q[2]);
						MINMAXH(q[3],q[1]);
						MINMAXH(q[3],q[2]);
						MINMAXH(p[0],q[0]);
						MINMAXH(p[1],q[1]);
						MINMAXH(p[2],q[2]);
						MINMAXH(p[3],q[3]);
					}//*/
				}
				ElemBlk[b][1]= j-1;
				++b; // compute number of blocks with compatible elements
				delete[]ElemConn;
            }
        }
        if( cfl>0.0 )
        {
            dt= cfl*hmin/Vmax();
            if( T==0.0 )
                T= dt*nt;
            else
            {
                nt= int(T/dt);
                if( dt*nt < T )
                    ++nt;
                T= dt*nt;
            }
        }
        // Display bounding rec and range for h
        if( DispInfo )
        {
            MPIPRINTF("\tBounding Box = (%g, %g, %g) - (%g, %g, %g)\n",
                      PT0.x, PT0.y, PT0.z, PT1.x, PT1.y, PT1.z);
            MPIPRINTF("\tMin(h) = %g, Max(h) = %g, dt = %g, nt = %i, T = %g, Vmax = %g\n",
                      hmin, hmax, dt, nt, T, Vmax());
        }

        // Close the file and clean up
        exoerr= ex_close(exoid);
        if( exoerr ) throw(3);

        delete[]Nodes_x;
        delete[]Nodes_y;
        delete[]Nodes_z;
        delete[]ElemBlkID;
    }
    catch(int e)
    {
        if( e==1 )
            {MPIPRINTF("\nError [%i] loading Exodus II mesh file [%s]\n", exoid, meshfile);}
        else if( e==2 )
        {
            if( exoerr<0 )
                {MPIPRINTF("\nError [%i] reading Exodus II mesh file [%s]\n", exoerr, meshfile);}
            else if(exoerr>0)
                MPIPRINTF("\nWarning [%i] reading Exodus II mesh file [%s]\n", exoerr, meshfile);
        }
        else if( e==3 )
        {
            if( exoerr<0 )
                {MPIPRINTF("\nError [%i] closing Exodus II mesh file [%s]\n", exoerr, meshfile);}
            else if(exoerr>0)
                MPIPRINTF("\nWarning [%i] closing Exodus II mesh file [%s]\n", exoerr, meshfile);
        }
        else if( e==4 )
        {
        	MPIPRINTF("\nError reading Exodus II file %s, no usable element blocks\n",meshfile);
        }
        else
        {
        	MPIPRINTF("\nError reading Exodus II file %s, (unclasified error e=%i)\n",meshfile, e);
        }
        return false;
    }
    return true;
#else
    MPIPRINTF("Error reading file %s, no Exodus II library\n", meshfile);
    return false;
#endif
}

//---------------------------------------------------------------------------
// Map the point (x0, z0) to a reference rectangle using a bilinera map
bool cDomain::InvMap(double x0, double z0, double &psi, double &eta, dVec3D *pt)
{
    const double T=0.01;
    try
    {
        double a,b,c,d,u,v;

        if( nDim!=2 )
            throw(1);
        u= pt[0].x-pt[1].x-pt[2].x+pt[3].x;
        v= pt[0].z-pt[1].z-pt[2].z+pt[3].z;
        a= u*(pt[0].z-pt[2].z) - v*(pt[0].x-pt[2].x);
        b= u*(z0-pt[0].z) - v*(x0-pt[0].x) - (pt[0].z-pt[2].z)*(pt[0].x-pt[1].x) + (pt[0].x-pt[2].x)*(pt[0].z-pt[1].z);
        c= (pt[0].z-pt[1].z)*(x0-pt[0].x) - (pt[0].x-pt[1].x)*(z0-pt[0].z);
        if( is_small(a) )
        {// linear mapping
            eta= -c/b;
            if( !InUnitInterval(eta,T) )
                throw(1); // the two solutions are out of bounds
        }
        else
        {// quadratic formula for eta
            d= sqr(b) - 4.0*a*c;
            if( d<0.0 )// complex solutions for eta
                throw(1);
            eta= ( sqrt(d) - b)/(2.0*a);
            if( !InUnitInterval(eta,T) )
            {// eta out of bounds, find the other solution
                eta= ( - sqrt(d) - b)/(2.0*a);
                if( !InUnitInterval(eta,T) )
                    throw(1); // the two solutions are out of bounds
            }
        }
        c= eta*u-pt[0].x+pt[1].x;
        if( is_small(c) )// check division by zero
            throw(1);
        psi= ( eta*(pt[0].x-pt[2].x) + x0 - pt[0].x)/c;
        if( InUnitInterval(psi,T) )
            return true;
    }
    catch(...)
    {
        psi= eta= -1.0;
        //MPIPRINTF("cDomain::InvMap: unable to map (%f, %f)\t",x0,z0);
    }
    return false;
}
//---------------------------------------------------------------------------
// Map the vector x to a reference rectangle using a trilinear mapping
bool cDomain::InvMap(dVec3D &x0, dVec3D &psi, dVec3D *pt)
{
    const double T=0.05;
    try
    {
        cAlgebra alg;
        dmat A(3,3,0.0), E(3,8,0.0);
        dvec X(3,0.0), Y(3,0.0);
        double det;
        if( nDim!=3 )
            throw(1);
        X[0]= x0.x - pt[0].x;
        X[1]= x0.y - pt[0].y;
        X[2]= x0.z - pt[0].z;
        // Copy the element nodes
        for( int i=0; i<3; ++i )
        {
            A[i][0]= pt[1][i] - pt[0][i];
            A[i][1]= pt[2][i] - pt[0][i];
            A[i][2]= pt[4][i] - pt[0][i];

            E[i][0]= pt[0][i];
            E[i][1]= pt[1][i] - pt[0][i];
            E[i][2]= pt[2][i] - pt[0][i];
            E[i][3]= pt[4][i] - pt[0][i];
            E[i][4]= pt[0][i] - pt[1][i] - pt[2][i] + pt[3][i];
            E[i][5]= pt[0][i] - pt[1][i] - pt[4][i] + pt[5][i];
            E[i][6]= pt[0][i] - pt[2][i] - pt[4][i] + pt[6][i];
        }
        if( alg.Gauss(A,X,Y)==false )
            throw(2);
        // First order aproximation
		psi.x = Y[0];
		psi.y = Y[1];
		psi.z = Y[2];
		// Second order correction
		det= TSP(1,2,3);
		if( is_small(det) )
			throw(2);
		psi.x-= (TSP(6,2,3)*Y[0]*Y[1] + TSP(5,2,3)*Y[0]*Y[2] + TSP(4,2,3)*Y[1]*Y[2])/det;
		psi.y-= (TSP(1,6,3)*Y[0]*Y[1] + TSP(1,5,3)*Y[0]*Y[2] + TSP(1,4,3)*Y[1]*Y[2])/det;
		psi.z-= (TSP(1,2,6)*Y[0]*Y[1] + TSP(1,2,5)*Y[0]*Y[2] + TSP(1,2,4)*Y[1]*Y[2])/det;
		// Check results
		if( !( InUnitInterval(psi.x,T) && InUnitInterval(psi.y,T) && InUnitInterval(psi.z,T) ) )
			throw(3);
    }
    catch(...)
    {
        psi.Init(-1.0);
        return false;
    }
    return true;
}

//---------------------------------------------------------------------------
// Map the point (x0, z0) to a reference rectangle using a bilinera map and a block
// a bilinear mapping
bool cDomain::InvMap(double x0, double z0, double &psi, double &eta, int b)
{
	if( InvMap(x0, z0, psi, eta, blk[b].pt) )
		return true;
	else
		MPIPRINTF("\ncDomain::InvMap: unable to map the point (%f, %f) to block %i\n",x0,z0, b);
	return false;
}
//---------------------------------------------------------------------------
// Map the vector x in element e to the master element using
// a bilinear mapping
bool cDomain::InvMap(dVec3D &x0, dVec3D &psi, int b)
{
	if( InvMap(x0, psi, blk[b].pt) )
		return true;
	else
		MPIPRINTF("\ncDomain::InvMap: unable to map the point (%f, %f, %f) to block %i\n", x0.x, x0.y, x0.z, b);
	return false;
}
//---------------------------------------------------------------------------
// Map the vector x in element e to the master element using
// a bilinear mapping
bool cDomain::InvMap(double x, double y, double z, double &psi, double &eta, double &zeta, int b)
{
    try
    {
        cAlgebra alg;
        dmat A(3,3,0.0), E(3,8,0.0);
        dvec X(3,0.0), Y(3,0.0);
        double det;
        if( nDim!=3 )
            throw(1);
        X[0]= x - blk[b].pt[0].x;
        X[1]= y - blk[b].pt[0].y;
        X[2]= z - blk[b].pt[0].z;
        // Copy the element nodes
        for( int i=0; i<3; ++i )
        {
            A[i][0]= blk[b].pt[1][i] - blk[b].pt[0][i];
            A[i][1]= blk[b].pt[2][i] - blk[b].pt[0][i];
            A[i][2]= blk[b].pt[4][i] - blk[b].pt[0][i];

            E[i][0]= blk[b].pt[0][i];
            E[i][1]= blk[b].pt[1][i] - blk[b].pt[0][i];
            E[i][2]= blk[b].pt[2][i] - blk[b].pt[0][i];
            E[i][3]= blk[b].pt[4][i] - blk[b].pt[0][i];
            E[i][4]= blk[b].pt[0][i] - blk[b].pt[1][i] - blk[b].pt[2][i] + blk[b].pt[3][i];
            E[i][5]= blk[b].pt[0][i] - blk[b].pt[1][i] - blk[b].pt[4][i] + blk[b].pt[5][i];
            E[i][6]= blk[b].pt[0][i] - blk[b].pt[2][i] - blk[b].pt[4][i] + blk[b].pt[6][i];
        }
        if( alg.Gauss(A,X,Y)==false )
            throw(2);
        // First order aproximation
		psi = Y[0];
		eta = Y[1];
		zeta= Y[2];
		// Second order correction
		det= TSP(1,2,3);
		if( is_small(det) )
			throw(2);
		psi -= (TSP(6,2,3)*Y[0]*Y[1] + TSP(5,2,3)*Y[0]*Y[2] + TSP(4,2,3)*Y[1]*Y[2])/det;
		eta -= (TSP(1,6,3)*Y[0]*Y[1] + TSP(1,5,3)*Y[0]*Y[2] + TSP(1,4,3)*Y[1]*Y[2])/det;
		zeta-= (TSP(1,2,6)*Y[0]*Y[1] + TSP(1,2,5)*Y[0]*Y[2] + TSP(1,2,4)*Y[1]*Y[2])/det;
		// Check results
		if( !( InUnitInterval(psi,0.05) && InUnitInterval(eta,0.05) && InUnitInterval(zeta,0.05) ) )
			throw(3);
    }
    catch(...)
    {
        psi= eta= zeta= -1.0;
        return false;
    }
    return true;
}

// Evaluate if a point is in a block
bool cDomain::InBlock(double x, double z, int b)
{
    double chi, eta;
    return InvMap(x,z,chi,eta,blk[b].pt);
}

// Evaluate if a point is in a block
bool cDomain::InBlock(double x, double y, double z, int b)
{
    dVec3D x0, psi;
	x0.Init(x,y,z);
    return InvMap(x0, psi, blk[b].pt);
}

//---------------------------------------------------------------------------
// Check if the point is inside the domain
bool cDomain::InDomain(double x, double z)
{
	double psi, eta;
	return InvMap(x, z, psi, eta, PT);//DD.VPT);
}
//---------------------------------------------------------------------------
// Check if the point is inside the domain
bool cDomain::InDomain(double x, double y, double z)
{
    dVec3D x0, psi;
	x0.Init(x,y,z);
	return InvMap(x0, psi, DD.VPT);
}

// Evaluate if a point is in a block
bool cDomain::InFracZone(double x, double z)
{
    double chi, eta;
    return InvMap(x,z,chi,eta,FZ);
}

// Evaluate if a point is in a block
bool cDomain::InFracZone(dVec3D &x)
{
    dVec3D psi;
    return InvMap(x, psi, FZ);
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::CalcMediaPar(int ix, int iz, double x, double z, double &vp0, double &vs0, double &rho0)
{
    for( int b=0; b<nBlocks; ++b )
    {
        if( InBlock(x,z,b) )
        //if(z<=blk[b].pt[3].z)
        {
            if( blk[b].Rho!=0.0 )
            {
                vp0 = blk[b].Vp;
                vs0 = blk[b].Vs;
                rho0= blk[b].Rho;
            }
            else// Non-homogeneous model
            {   // Check indices
                if( ix<0 ) ix= 0;
                if( iz<0 ) iz= 0;
                if( ix>=n.x ) ix= n.x-1;
                if( iz>=n.z ) iz= n.z-1;
                // Get parameters
                vp0 = Vp [iz][ix];
                vs0 = Vs [iz][ix];
                rho0= Rho[iz][ix];
            }
            return true;
        }
    }
    return false;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::CalcMediaPar(double x, double z, double &vp0, double &vs0, double &rho0)
{
    int b;
    for( b=0; b<nBlocks; ++b )
    {
        if( InBlock(x,z,b) )
        //if(z<=blk[b].pt[3].z)
        {
            CalcMediaPar(x,z,vp0,vs0,rho0,b);
            break;
        }
    }
    if( b>=nBlocks )
    {
        MPIPRINTF("\ncDomain::CalcMediaPar: unable to find the block for (%g,%g)",x,z);
        return false;
    }
    return true;
}

// Interpolate the media parameters at the point (x,y,z)
bool cDomain::CalcMediaPar(double x, double y, double z,
        double &vp0, double &vs0, double &rho0)
{
    int b;
    for( b=0; b<nBlocks; ++b )
    {
        if( InBlock(x,y,z,b) )
        //if(z<=blk[b].pt[7].z)
        {
            CalcMediaPar(x,y,z,vp0,vs0,rho0,b);
            break;
        }
    }
    if( b>=nBlocks )
    {
        MPIPRINTF("\ncDomain::CalcMediaPar: unable to find the block for (%g,%g,%g)",
                x,y,z);
        return false;
    }
    return true;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::CalcMediaPar(double x, double z, double &vp0, double &vs0, double &ro, int b)
{
    int n1, n2;

    if( nBlocks<0 || b<0 )
        return false;

    if( blk[b].Rho!=0.0 )
    {
        vp0 = blk[b].Vp;
        vs0 = blk[b].Vs;
        ro  = blk[b].Rho;
    }
    else
    {   // Non-homogeneous model
        if( Vp.IsEmpty() || Vs.IsEmpty() || Rho.IsEmpty() )
            return false;
        if( x>PT1.x || x<PT0.x || z>PT1.z || z<PT0.z )
            return false;

        n1 = int(floor(n.x*(x - PT0.x)/(PT1.x - PT0.x)));
        n2 = int(floor(n.z*(z - PT0.z)/(PT1.z - PT0.z)));
        if( n1 >= n.x ) n1= n.x-1;
        if( n2 >= n.z ) n2= n.z-1;

        vp0 = Vp[n2][n1];
        vs0 = Vs[n2][n1];
        ro  = Rho[n2][n1];
    }
    return true;
}

// Interpolate the media parameters at the point (x,y,z)
bool cDomain::CalcMediaPar(double x, double y, double z,
        double &vp0, double &vs0, double &ro, int b)
{
    if( nBlocks<0 || b<0 )
        return false;

    if( blk[b].Rho<=0.0 )
        return false;

    vp0 = blk[b].Vp;
    vs0 = blk[b].Vs;
    ro  = blk[b].Rho;

    return true;
}

// Interpolate the media parameters at the point (x,y,z)
bool cDomain::CalcMediaPar(double &vp0, double &vs0, double &ro, int b)
{
    if( nBlocks<0 || b<0 )
        return false;

    if( blk[b].Rho<=0.0 )
        return false;

    vp0 = blk[b].Vp;
    vs0 = blk[b].Vs;
    ro  = blk[b].Rho;

    return true;
}

// Get the media parameters at the point (x,z)
bool cDomain::GetRho(double &rho, int b)
{
    if( nBlocks<0 || b<0 || blk[b].Rho==0.0 )
        return false;

    rho= blk[b].Rho;
    return true;
}

// Get the media parameters at the point (x,z)
bool cDomain::CalcMediaLame(int ix, int iz, double x, double z, double &la, double &mu, double &r)
{
    double vp0, vs0;

    if( CalcMediaPar(ix,iz, x, z,vp0,vs0,r) )
    {
        mu= r*sqr(vs0);
        la= r*sqr(vp0) - 2.0*mu;
        return true;
    }
    return false;
}

// Get the media parameters at the point (x,y,z)
bool cDomain::CalcMediaLame(int ix, int iy, int iz,
        double &la, double &mu, double &r)
{
    double vp0, vs0;

    if( CalcMediaPar(ix,iy,iz,vp0,vs0,r) )
    {
        mu= r*sqr(vs0);
        la= r*sqr(vp0) - 2.0*mu;
        return true;
    }
    return false;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::CalcMediaLame(double x, double z, double &la, double &mu, double &r)
{
    double vp0, vs0;

    if( CalcMediaPar(x,z,vp0,vs0,r) )
    {
        mu= r*sqr(vs0);
        la= r*sqr(vp0) - 2.0*mu;
        return true;
    }
    return false;
}

// Interpolate the media parameters at the point (x,y,z)
bool cDomain::CalcMediaLame(double x, double y, double z,
        double &la, double &mu, double &r)
{
    double vp0, vs0;

    if( CalcMediaPar(x,y,z,vp0,vs0,r) )
    {
        mu= r*sqr(vs0);
        la= r*sqr(vp0) - 2.0*mu;
        return true;
    }
    return false;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::CalcMediaLame(double x, double z, double &la, double &mu, double &r, int b)
{
    double vp0, vs0;

    if( CalcMediaPar(x,z,vp0,vs0,r,b) )
    {
        mu= r*sqr(vs0);
        la= r*sqr(vp0) - 2.0*mu;
        return true;
    }
    return false;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::CalcMediaLame(double x, double y, double z,
        double &la, double &mu, double &r, int b)
{
    double vp0, vs0;

    if( CalcMediaPar(x,y,z,vp0,vs0,r,b) )
    {
        mu= r*sqr(vs0);
        la= r*sqr(vp0) - 2.0*mu;
        return true;
    }
    return false;
}

// Interpolate the media parameters of block b
bool cDomain::CalcMediaLame(double &la, double &mu, double &r, int b)
{
    double vp0, vs0;

    if( CalcMediaPar(vp0,vs0,r,b) )
    {
        mu= r*sqr(vs0);
        la= r*sqr(vp0) - 2.0*mu;
        return true;
    }
    return false;
}

// Get the elastic tensor of block b
cElasTensor *cDomain::CalcMediaCij(double &r, int b)
{
    if( blk[b].C!=NULL )
    {
    	r= blk[b].Rho;
        return blk[b].C;
    }
    return NULL;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::FwdMap(double psi, double eta, double &x, double &z)
{
    double mapx[2], mapz[2], map[4];

    mapx[0]= 1.0 - psi;
    mapx[1]= psi;
    mapz[0]= 1.0 - eta;
    mapz[1]= eta;

    map[0]= mapx[0]*mapz[0];
    map[1]= mapx[1]*mapz[0];
    map[2]= mapx[0]*mapz[1];
    map[3]= mapx[1]*mapz[1];

    x= z= 0.0;
    for( int i=0; i<4; ++i )
    {
        x+= map[i]*PT[i].x;
        z+= map[i]*PT[i].z;
    }
    return true;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::FwdMap(dVec3D &psi, dVec3D &x)
{
    double mapx[2], mapy[2], mapz[2], map[8];

    mapx[0]= 1.0 - psi.x;
    mapx[1]= psi.x;
    mapy[0]= 1.0 - psi.y;
    mapy[1]= psi.y;
    mapz[0]= 1.0 - psi.z;
    mapz[1]= psi.z;

    map[0]= mapx[0]*mapy[0]*mapz[0];
    map[1]= mapx[1]*mapy[0]*mapz[0];
    map[2]= mapx[0]*mapy[1]*mapz[0];
    map[3]= mapx[1]*mapy[1]*mapz[0];
    map[4]= mapx[0]*mapy[0]*mapz[1];
    map[5]= mapx[1]*mapy[0]*mapz[1];
    map[6]= mapx[0]*mapy[1]*mapz[1];
    map[7]= mapx[1]*mapy[1]*mapz[1];

    x.Init(0.0);
    for( int i=0; i<4; ++i )
    {
        x+= PT[i]*map[i];
    }
    return true;
}

// Inveres bilinera map from the domain to a reference square
bool cDomain::InvMap(double x, double z, double &psi, double &eta)
{
	return InvMap(x, z, psi, eta, PT);
}

// Inveres bilinera map from the domain to a reference square
bool cDomain::InvMap(dVec3D &x0, dVec3D &psi)
{
	return InvMap(x0, psi, PT);
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::GetLamAcoust(double x, double z, double &lambda, double &r, int b)
{
    double aux;

    if( CalcMediaPar(x,z,lambda,aux,r,b) )
    {
        lambda= r*sqr(lambda);
        return true;
    }
    return false;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::GetLamAcoust(double x, double y, double z, double &lambda, double &r, int b)
{
    double aux;

    if( CalcMediaPar(x,y,z,lambda,aux,r,b) )
    {
        lambda= r*sqr(lambda);
        return true;
    }
    return false;
}

// Interpolate the media parameters at the point (x,z)
bool cDomain::GetLamAcoust(double &lambda, double &r, int b)
{
    double aux;

    if( CalcMediaPar(lambda,aux,r,b) )
    {
        lambda= r*sqr(lambda);
        return true;
    }
    return false;
}

// Resample the media parameters in an n.x by n.z mesh
void cDomain::Resample(int mx, int mz)
{
    cMat2<float> newVp(mx,mz,0.0), newVs(mx,mz,0.0), newRho(mx,mz,0.0);
    double hx, hz, vp0, vs0, r0;

    // Return if the requested is the same as the actual
    if( mx==n.x && mz==n.z )
        return;

    hx= (PT1.x - PT0.x)/(mx - 1);
    hz= (PT1.z - PT0.z)/(mz - 1);

    for( int i=0; i<mz; ++i )
        for( int j=0; j<mx; ++j )
        {
            CalcMediaPar(hx*j, hz*i, vp0, vs0, r0);
            newVp[i][j] = (float) vp0;
            newVs[i][j] = (float) vs0;
            newRho[i][j]= (float) r0;
        }
    Vp = newVp;
    Vs = newVs;
    Rho= newRho;
    n.x= mx;
    n.z= mz;
    //del.x= (float) hx;
    //del.z= (float) hz;
}

// Preprocessor for the Domain Decomposition
// TO-DO: reimplement the domain decomposition for a general quadrilateral
void cDomain::Preproc(void)
{
#ifdef ENABLEMETIS
    try
	{
        MPI_id= MPI::COMM_WORLD.Get_rank();
        MPI_np= MPI::COMM_WORLD.Get_size();

        DD.PT0= DD.VPT0= PT0;
        DD.PT1= DD.VPT1= PT1;

        for( int i=0; i<8; ++i )
            DD.VPT[i]= DD.PT[i]= PT[i];

        // Setup the resolution of the snapshots
        GlobalRes.x= Res.x;
        GlobalRes.y= int(GlobalRes.x*(PT1.y - PT0.y)/(PT1.x - PT0.x));
        GlobalRes.z= int(GlobalRes.x*(PT1.z - PT0.z)/(PT1.x - PT0.x));

        Res.x= int(GlobalRes.x*(DD.VPT1.x - DD.VPT0.x)/(PT1.x - PT0.x));
        Res.y= int(GlobalRes.x*(DD.VPT1.y - DD.VPT0.y)/(PT1.x - PT0.x));
        Res.z= int(GlobalRes.x*(DD.VPT1.z - DD.VPT0.z)/(PT1.x - PT0.x));

        Org.x= int(GlobalRes.x*(PT0.x + DD.VPT0.x)/(PT1.x - PT0.x));
        Org.y= int(GlobalRes.x*(PT0.y + DD.VPT0.y)/(PT1.x - PT0.x));
        Org.z= int(GlobalRes.x*(PT0.z + DD.VPT0.z)/(PT1.x - PT0.x));
	}
    catch(int err)
    {
        printf("MPI Error %i\n",err);
    }

#elif defined ENABLEMPI
    try
    {
        double delta= (PT1.x - PT0.x)/nElem.x;

        MPI_id= MPI::COMM_WORLD.Get_rank();
        MPI_np= MPI::COMM_WORLD.Get_size();

        if(MPI_np<=1)
        {
            DD.PT0= DD.VPT0= PT0;
            DD.PT1= DD.VPT1= PT1;

            for( int i=0; i<8; ++i )
                DD.VPT[i]= DD.PT[i]= PT[i];
        }
        else
        {
            int DD_nElem0= nElem.x/MPI_np;
            DD.VPT0.x= PT0.x + delta*DD_nElem0*MPI_id;

            if( MPI_id==0 )
            {
                DD.PT0.x= DD.VPT0.x;
                DD.VPT1.x= DD.VPT0.x + delta*DD_nElem0;
                ++DD_nElem0;
            }
            else if( MPI_id == MPI_np-1 )
            {
                DD.VPT1.x= PT1.x;
                DD_nElem0= nElem.x - DD_nElem0*MPI_id + 1;
                DD.PT0.x= PT1.x - delta*DD_nElem0;
            }
            else
            {
                DD.PT0.x= DD.VPT0.x - delta;
                DD.VPT1.x= DD.VPT0.x + delta*DD_nElem0;
                DD_nElem0+=2;
            }
            DD.PT1.x= DD.PT0.x + delta*DD_nElem0;
            DD.PT0.y= DD.VPT0.y= PT0.y;
            DD.PT0.z= DD.VPT0.z= PT0.z;
            DD.PT1.y= DD.VPT1.y= PT1.y;
            DD.PT1.z= DD.VPT1.z= PT1.z;
            nElem.x= DD_nElem0;
            if( tRandFrac && !bFracZone)
            {
                int f= nFracs/MPI_np;
                if( MPI_id==0 )
                    nFracs= nFracs - f*(MPI_np-1);
                else
                    nFracs= f;
            }
        }
        nElem.z= 0;
        for(int b=0; b<nBlocks; ++b )
        {
            double m1, m2, a;
            if( blk[b].Type!=BT_LAYER )
                continue;

            // Number of elements in each block
            blk[b].nElem.x= nElem.x;
            nElem.z+= blk[b].nElem.z;
            if( nDim==2 )
            {// Adjust the block, top
                m1= (blk[b].pt[1].z - blk[b].pt[0].z)/(blk[b].pt[1].x - blk[b].pt[0].x);
                a= blk[b].pt[0].z - m1*blk[b].pt[0].x;
                blk[b].pt[0].x= DD.PT0.x;
                blk[b].pt[1].x= DD.PT1.x;
                blk[b].pt[0].z= m1*DD.PT0.x + a;
                blk[b].pt[1].z= m1*DD.PT1.x + a;
                // Adjust the block, bottom
                m1= (blk[b].pt[3].z - blk[b].pt[2].z)/(blk[b].pt[3].x - blk[b].pt[2].x);
                a= blk[b].pt[2].z - m1*blk[b].pt[2].x;
                blk[b].pt[2].x= DD.PT0.x;
                blk[b].pt[3].x= DD.PT1.x;
                blk[b].pt[2].z= m1*DD.PT0.x + a;
                blk[b].pt[3].z= m1*DD.PT1.x + a;
            }
            else // 3D
            {// Adjust the blocks, top
                a=    blk[b].pt[0].y*(blk[b].pt[2].x - blk[b].pt[1].x)
                    + blk[b].pt[1].y*(blk[b].pt[0].x - blk[b].pt[2].x)
                    + blk[b].pt[2].y*(blk[b].pt[1].x - blk[b].pt[0].x);

                m1= ( blk[b].pt[0].y*(blk[b].pt[2].z - blk[b].pt[1].z)
                    + blk[b].pt[1].y*(blk[b].pt[0].z - blk[b].pt[2].z)
                    + blk[b].pt[2].y*(blk[b].pt[1].z - blk[b].pt[0].z) )/a;

                m2= ( blk[b].pt[0].z*(blk[b].pt[2].x - blk[b].pt[1].x)
                    + blk[b].pt[1].z*(blk[b].pt[0].x - blk[b].pt[2].x)
                    + blk[b].pt[2].z*(blk[b].pt[1].x - blk[b].pt[0].x) )/a;

                a= blk[b].pt[0].z - m1*blk[b].pt[0].x - m2*blk[b].pt[0].y;

                blk[b].pt[0].x= DD.PT0.x;
                blk[b].pt[1].x= DD.PT1.x;
                blk[b].pt[2].x= DD.PT0.x;
                blk[b].pt[3].x= DD.PT1.x;

                blk[b].pt[0].y= DD.PT0.y;
                blk[b].pt[1].y= DD.PT0.y;
                blk[b].pt[2].y= DD.PT1.y;
                blk[b].pt[3].y= DD.PT1.y;

                blk[b].pt[0].z= m1*DD.PT0.x + m2*DD.PT0.y + a;
                blk[b].pt[1].z= m1*DD.PT1.x + m2*DD.PT0.y + a;
                blk[b].pt[2].z= m1*DD.PT0.x + m2*DD.PT1.y + a;
                blk[b].pt[3].z= m1*DD.PT1.x + m2*DD.PT1.y + a;

                a=    blk[b].pt[4].y*(blk[b].pt[6].x - blk[b].pt[5].x)
                    + blk[b].pt[5].y*(blk[b].pt[4].x - blk[b].pt[6].x)
                    + blk[b].pt[6].y*(blk[b].pt[5].x - blk[b].pt[4].x);

                m1= ( blk[b].pt[4].y*(blk[b].pt[6].z - blk[b].pt[5].z)
                    + blk[b].pt[5].y*(blk[b].pt[4].z - blk[b].pt[6].z)
                    + blk[b].pt[6].y*(blk[b].pt[5].z - blk[b].pt[4].z) )/a;

                m2= ( blk[b].pt[4].z*(blk[b].pt[6].x - blk[b].pt[5].x)
                    + blk[b].pt[5].z*(blk[b].pt[4].x - blk[b].pt[6].x)
                    + blk[b].pt[6].z*(blk[b].pt[5].x - blk[b].pt[4].x) )/a;

                a= blk[b].pt[4].z - m1*blk[b].pt[4].x - m2*blk[b].pt[4].y;

                blk[b].pt[4].x= DD.PT0.x;
                blk[b].pt[5].x= DD.PT1.x;
                blk[b].pt[6].x= DD.PT0.x;
                blk[b].pt[7].x= DD.PT1.x;

                blk[b].pt[4].y= DD.PT0.y;
                blk[b].pt[5].y= DD.PT0.y;
                blk[b].pt[6].y= DD.PT1.y;
                blk[b].pt[7].y= DD.PT1.y;

                blk[b].pt[4].z= m1*DD.PT0.x + m2*DD.PT0.y + a;
                blk[b].pt[5].z= m1*DD.PT1.x + m2*DD.PT0.y + a;
                blk[b].pt[6].z= m1*DD.PT0.x + m2*DD.PT1.y + a;
                blk[b].pt[7].z= m1*DD.PT1.x + m2*DD.PT1.y + a;
            }
        }
        if( nDim==2 )
            nElem.y= nElem.z;

        // Setup the resolution of the snapshots
        GlobalRes.x= Res.x;
        GlobalRes.y= int(GlobalRes.x*(PT1.y - PT0.y)/(PT1.x - PT0.x));
        GlobalRes.z= int(GlobalRes.x*(PT1.z - PT0.z)/(PT1.x - PT0.x));

        Res.x= int(GlobalRes.x*(DD.VPT1.x - DD.VPT0.x)/(PT1.x - PT0.x));
        Res.y= int(GlobalRes.x*(DD.VPT1.y - DD.VPT0.y)/(PT1.x - PT0.x));
        Res.z= int(GlobalRes.x*(DD.VPT1.z - DD.VPT0.z)/(PT1.x - PT0.x));

        Org.x= int(GlobalRes.x*(-PT0.x + DD.VPT0.x)/(PT1.x - PT0.x));
        Org.y= int(GlobalRes.x*(-PT0.y + DD.VPT0.y)/(PT1.x - PT0.x));
        Org.z= int(GlobalRes.x*(-PT0.z + DD.VPT0.z)/(PT1.x - PT0.x));

        if( MPI_id==0 )
        {
            MPI_Res.Init(MPI_np);
            MPI_Org.Init(MPI_np);
        }
        MPI::COMM_WORLD.Gather(&Res, 3, MPI::INT, MPI_Res.Data(), 3, MPI::INT, 0);
        MPI::COMM_WORLD.Gather(&Org, 3, MPI::INT, MPI_Org.Data(), 3, MPI::INT, 0);

        // Test code
        dVec3D *buff_pt0= new dVec3D[MPI_np],
        *buff_pt1= new dVec3D[MPI_np];

        MPI::COMM_WORLD.Gather(&DD.PT0, 3, MPI::DOUBLE, buff_pt0, 3, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Gather(&DD.PT1, 3, MPI::DOUBLE, buff_pt1, 3, MPI::DOUBLE, 0);
        /*
        printf("[P%.2i] %i %i %i\n", MPI_id, Res.x, Res.y, Res.z);
        if( MPI_id==0 )
            for(int i=0; i<MPI_np; ++i)
            {
                printf("[P%.2i] Res [%i %i %i] ", i, MPI_Res[i].x, MPI_Res[i].y, MPI_Res[i].z);
                printf("Org [%i %i %i] ", MPI_Org[i].x, MPI_Org[i].y, MPI_Org[i].z);
                printf("pt0 [%g %g %g] - ", buff_pt0[i].x, buff_pt0[i].y, buff_pt0[i].z);
                printf("pt1 [%g %g %g]\n", buff_pt1[i].x, buff_pt1[i].y, buff_pt1[i].z);
            }
        */
        delete[]buff_pt0;
        delete[]buff_pt1;
        //*/
        //MPI::COMM_WORLD.Barrier();
    }
    catch(int err)
    {
        printf("MPI Error %i\n",err);
    }
#else
    DD.PT0.x= DD.VPT0.x= PT0.x;
    DD.PT0.y= DD.VPT0.y= PT0.y;
    DD.PT0.z= DD.VPT0.z= PT0.z;
    DD.PT1.x= DD.VPT1.x= PT1.x;
    DD.PT1.y= DD.VPT1.y= PT1.y;
    DD.PT1.z= DD.VPT1.z= PT1.z;

    for( int i=0; i<8; ++i )
        DD.VPT[i]= DD.PT[i]= PT[i];

    // Setup the resolution of the snapshots
    GlobalRes.x= Res.x;
    GlobalRes.y= int(GlobalRes.x*(PT1.y - PT0.y)/(PT1.x - PT0.x));
    GlobalRes.z= int(GlobalRes.x*(PT1.z - PT0.z)/(PT1.x - PT0.x));

    Res.x= int(GlobalRes.x*(DD.VPT1.x - DD.VPT0.x)/(PT1.x - PT0.x));
    Res.y= int(GlobalRes.x*(DD.VPT1.y - DD.VPT0.y)/(PT1.x - PT0.x));
    Res.z= int(GlobalRes.x*(DD.VPT1.z - DD.VPT0.z)/(PT1.x - PT0.x));

    Org.x= int(GlobalRes.x*(PT0.x + DD.VPT0.x)/(PT1.x - PT0.x));
    Org.y= int(GlobalRes.x*(PT0.y + DD.VPT0.y)/(PT1.x - PT0.x));
    Org.z= int(GlobalRes.x*(PT0.z + DD.VPT0.z)/(PT1.x - PT0.x));
#endif
}
