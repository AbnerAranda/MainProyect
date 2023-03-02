/*
 * Declaration of the domain class to manage the physical domain,
 * read the input files and store all the input parameters.
 * Jonas D. De Basabe (jonas@cicese.mx)
 * File:   cdomain.{cpp,h}
 */
//---------------------------------------------------------------------------
#ifndef CDOMAINHDR
#define CDOMAINHDR

#include "macros.h"
#include "celastensor.h"

class cEleMesh;
//---------------------------------------------------------------------------
class cDomain
{
protected:
    bool CheckBlocks(void);
    void ExtractParam(char *buff);
    void InitBlocks(void);
    bool ReadMediaPar(char *szfile, cMat2<float> &a);
    bool ReadMediaParBin(char *szfile, cMat2<float> &a);
    bool ReadMediaParTxt(char *szfile, cMat2<float> &a);
    bool ReadMeshTxt(dmat &Nod, cEleMesh &Ele, cMat2<int> &ElemBlk);
    bool ReadMeshBin(dmat &Nod, cEleMesh &Ele, cMat2<int> &ElemBlk);
    bool InvMap(double x0, double z0, double &psi, double &eta, dVec3D *pt);
    bool InvMap(double x0, double y0, double z0, double &psi, double &eta, double &zeta,  double *x, double *y, double *z);
    bool InvMap(dVec3D &x0, dVec3D &psi, dVec3D *pt);
    bool InvMap(double x, double z, double &psi, double &eta, int b);
    bool InvMap(dVec3D &x0, dVec3D &psi, int b);
    bool InvMap(double x, double y, double z, double &psi, double &eta, double &zeta, int b);
    int ch2int(char c){return c - '0';}
    bool str2int(char *c, int &i)
    {
    	int e= sscanf(c,"%i",&i);
    	return (e>=1);
    }

public:
    struct DGparamters
    {// DG parameters
        bool bExactInt, bPenVel;
        double  Penalty;
        TDGMethod Type;
    } DG;

    struct EGparameters
    {// EG parameters
        int Ord;
        TBasis Type;
    } EG;

    struct Source
    {// Source parameters
        bool bSave;
        int nShots, tFunc;
        double  PkFq, tWidth, Amp, tensr[6];
        dVec3D loc, nor, sWidth;
        TSrcType Type;
    } src;

    struct Receivers
    {// Receivers parameters
        int nTraces, nSeismo;
        dVec3D r0, d;
        cVector<dVec3D> Seismo;
    } rec;

    struct subdomain
    {
        bool bHomog;
        iVec3D nElem;
        double Vp, Vs, Rho, Vmax, Penalty, FracZt, FracZn;
        cElasTensor *C;
        dVec3D pt[8];
        dVec3D h;
        int nFrac, tFrac;
        TBlockType Type;
    };

    struct DomainDecomp
    {
        dVec3D PT[8],   // Whole domain after decomposition
                VPT[8],   // Domain without overlaps
                PT0, PT1, 	// Whole domain after decomposition
                VPT0, VPT1; // Domain without overlaps
    } DD;

    struct fracture
    {
        dVec3D pt[4];   // fracture coordinates
        double Zt, Zn;  // fracture compliance
    };

    struct GrDisp
    {
        bool ContTime;  // Continuous time analysis
        bool Plot;      // Create the GNU plot
        int Type;       // Type of plot
    } gd;

    cVector<fracture> frac;
    cVector<subdomain> blk;

    bool bElast, bFreeSurf, bHomog, bSaveMesh, bFracZone, bSeismoDXU;
    int DispInfo, nBlocks, nFracs, nt,
        nTaper, ElemOrd, nTSave, nDim, MeshFac,
        MPI_id, MPI_np, RandSeed;
    double dt, T, cfl,
        tapalpha, tapbeta,
        vphom, vshom, rhohom, fTSave;// Vmax,
    dVec3D PT0, PT1, delfrac;
    dVec3D PT[8], FZ[8]; // points to define the whole domain and frac zone
    iVec3D Res, Org, GlobalRes, n, nElem;
    cVector<iVec3D> MPI_Org, MPI_Res;

    char *vpfile, *vsfile, *rhofile,
        *meshfile, *logfile, *outdir;

    cMat2<float> Rho, Vp, Vs;

    TExeMode ExeMode;
    TBasis BasisType;
    TModel Model;
    TMethod Method;
    TMeshType MeshType;
    TBoundary Boundary;
    TTSMethod TSMethod;
    TSnapshots Snapsh;
    TRandFract tRandFrac;

    cDomain(char *parfile=NULL); // cDomain() == cDomain(NULL)
    ~cDomain(void);

    void Init(void);
    bool PhysRec(void);
    bool InDomain(double x, double z);
    bool InDomain(double x, double y, double z);
    bool InFracZone(double x, double z);
    bool InFracZone(dVec3D &x);
    bool InBlock(double x, double z, int b);
    bool InBlock(double x, double y, double z, int b);
    bool GetXBlocks(cVector<int> &XBlks);
    double GetTimeInc(double h, double v);
    double GetMaxV(void);
    double GetMaxV(int b);
    void GetDelta(double &h, int &b);
    void GetDelta(double &h, double &v);
    int GetNTSave(void);
    bool SetFracEqMedia(void);
    void Display(void);
    bool ReadParFile(char *parfile);
    bool ReadMesh(cEleMesh &Ele, cMat2<int> &ElemBlk);

    // Get the media parameters by block number
    bool GetRho(double &rho, int b);
    bool CalcMediaPar(double &vp0, double &vs0, double &rho0, int b);
    bool CalcMediaLame(double &la, double &mu, double &ro, int b);
    bool GetLamAcoust(double &la, double &r, int b);
    // Get the media parameters by x, y and z index
    bool CalcMediaPar(int ix, int iz, double x, double z,
            double &vp0, double &vs0, double &rho0);
    bool CalcMediaLame(int ix, int iz, double x, double z,
            double &la, double &mu, double &ro);
    bool CalcMediaLame(int ix, int iy, int iz,
            double &la, double &mu, double &ro);
    // Get media parameters by coordinates
    bool CalcMediaPar(double x, double z,
            double &vp0, double &vs0, double &rho0);
    bool CalcMediaPar(double x, double y, double z,
            double &vp0, double &vs0, double &rho0);
    bool CalcMediaPar(double x, double z,
            double &vp0, double &vs0, double &rho0, int b);
    bool CalcMediaPar(double x, double y, double z,
            double &vp0, double &vs0, double &rho0, int b);
    bool CalcMediaLame(double x, double z,
            double &la, double &mu, double &ro);
    bool CalcMediaLame(double x, double y, double z,
            double &la, double &mu, double &ro);
    bool CalcMediaLame(double x, double z,
            double &la, double &mu, double &ro, int b);
    bool CalcMediaLame(double x, double y, double z,
            double &la, double &mu, double &ro, int b);
    cElasTensor *CalcMediaCij(double &r, int b);
    bool GetLamAcoust(double x, double z, double &lambda, double &r, int b);
    bool GetLamAcoust(double x, double y, double z, double &lambda, double &r, int b);
    // Forward and inverse mapping from a reference to the physical point using the domain quadrilateral
    bool FwdMap(double psi, double eta, double &x, double &z);
    bool InvMap(double x, double z, double &psi, double &eta);
    bool FwdMap(dVec3D &psi, dVec3D &x);
    bool InvMap(dVec3D &x, dVec3D &psi);
    //
    void Resample(int mx, int mz);
    void Preproc(void);
    //
    double Vmax(void)
    {// this does not take into account the models introduced by raw binary files
    	double v=0.0;
    	for(int i=0; i<nBlocks; ++i)
    	{
    		if(v<blk[i].Vp)
    			v= blk[i].Vp;
    	}
    	return v;
    }
};
#endif
