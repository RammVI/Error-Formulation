/* ./src_f77/dsyevx.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <punc/vf2c.h>

/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
	integer *iwork, integer *ifail, integer *info, ftnlen jobz_len, 
	ftnlen range_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nb, jj;
    static doublereal eps, vll, vuu, tmp1;
    static integer indd, inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static integer lopt, itmp1, indee;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lower, wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, indibl;
    static logical valeig;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal abstll, bignum;
    static integer indtau, indisp;
    extern /* Subroutine */ int dstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dsterf_(integer *, doublereal *, doublereal *, integer *);
    static integer indiwo, indwkn;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int dorgtr_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), dsteqr_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dormtr_(char *, char *, char *, integer *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer llwrkn, llwork, nsplit;
    static doublereal smlnum;
    extern /* Subroutine */ int dsytrd_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSYEVX computes selected eigenvalues and, optionally, eigenvectors */
/*  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be */
/*  selected by specifying either a range of values or a range of indices */
/*  for the desired eigenvalues. */

/*  Arguments */
/*  ========= */

/*  JOBZ    (input) CHARACTER*1 */
/*          = 'N':  Compute eigenvalues only; */
/*          = 'V':  Compute eigenvalues and eigenvectors. */

/*  RANGE   (input) CHARACTER*1 */
/*          = 'A': all eigenvalues will be found. */
/*          = 'V': all eigenvalues in the half-open interval (VL,VU] */
/*                 will be found. */
/*          = 'I': the IL-th through IU-th eigenvalues will be found. */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  Upper triangle of A is stored; */
/*          = 'L':  Lower triangle of A is stored. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N) */
/*          On entry, the symmetric matrix A.  If UPLO = 'U', the */
/*          leading N-by-N upper triangular part of A contains the */
/*          upper triangular part of the matrix A.  If UPLO = 'L', */
/*          the leading N-by-N lower triangular part of A contains */
/*          the lower triangular part of the matrix A. */
/*          On exit, the lower triangle (if UPLO='L') or the upper */
/*          triangle (if UPLO='U') of A, including the diagonal, is */
/*          destroyed. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  VL      (input) DOUBLE PRECISION */
/*  VU      (input) DOUBLE PRECISION */
/*          If RANGE='V', the lower and upper bounds of the interval to */
/*          be searched for eigenvalues. VL < VU. */
/*          Not referenced if RANGE = 'A' or 'I'. */

/*  IL      (input) INTEGER */
/*  IU      (input) INTEGER */
/*          If RANGE='I', the indices (in ascending order) of the */
/*          smallest and largest eigenvalues to be returned. */
/*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/*          Not referenced if RANGE = 'A' or 'V'. */

/*  ABSTOL  (input) DOUBLE PRECISION */
/*          The absolute error tolerance for the eigenvalues. */
/*          An approximate eigenvalue is accepted as converged */
/*          when it is determined to lie in an interval [a,b] */
/*          of width less than or equal to */

/*                  ABSTOL + EPS *   max( |a|,|b| ) , */

/*          where EPS is the machine precision.  If ABSTOL is less than */
/*          or equal to zero, then  EPS*|T|  will be used in its place, */
/*          where |T| is the 1-norm of the tridiagonal matrix obtained */
/*          by reducing A to tridiagonal form. */

/*          Eigenvalues will be computed most accurately when ABSTOL is */
/*          set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
/*          If this routine returns with INFO>0, indicating that some */
/*          eigenvectors did not converge, try setting ABSTOL to */
/*          2*DLAMCH('S'). */

/*          See "Computing Small Singular Values of Bidiagonal Matrices */
/*          with Guaranteed High Relative Accuracy," by Demmel and */
/*          Kahan, LAPACK Working Note #3. */

/*  M       (output) INTEGER */
/*          The total number of eigenvalues found.  0 <= M <= N. */
/*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */

/*  W       (output) DOUBLE PRECISION array, dimension (N) */
/*          On normal exit, the first M elements contain the selected */
/*          eigenvalues in ascending order. */

/*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M)) */
/*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/*          contain the orthonormal eigenvectors of the matrix A */
/*          corresponding to the selected eigenvalues, with the i-th */
/*          column of Z holding the eigenvector associated with W(i). */
/*          If an eigenvector fails to converge, then that column of Z */
/*          contains the latest approximation to the eigenvector, and the */
/*          index of the eigenvector is returned in IFAIL. */
/*          If JOBZ = 'N', then Z is not referenced. */
/*          Note: the user must ensure that at least max(1,M) columns are */
/*          supplied in the array Z; if RANGE = 'V', the exact value of M */
/*          is not known in advance and an upper bound must be used. */

/*  LDZ     (input) INTEGER */
/*          The leading dimension of the array Z.  LDZ >= 1, and if */
/*          JOBZ = 'V', LDZ >= max(1,N). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The length of the array WORK.  LWORK >= max(1,8*N). */
/*          For optimal efficiency, LWORK >= (NB+3)*N, */
/*          where NB is the max of the blocksize for DSYTRD and DORMTR */
/*          returned by ILAENV. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  IWORK   (workspace) INTEGER array, dimension (5*N) */

/*  IFAIL   (output) INTEGER array, dimension (N) */
/*          If JOBZ = 'V', then if INFO = 0, the first M elements of */
/*          IFAIL are zero.  If INFO > 0, then IFAIL contains the */
/*          indices of the eigenvectors that failed to converge. */
/*          If JOBZ = 'N', then IFAIL is not referenced. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, then i eigenvectors failed to converge. */
/*                Their indices are stored in array IFAIL. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;

    /* Function Body */
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;

    *info = 0;
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else {
	if (valeig) {
	    if (*n > 0 && *vu <= *vl) {
		*info = -8;
	    }
	} else if (indeig) {
	    if (*il < 1 || *il > max(1,*n)) {
		*info = -9;
	    } else if (*iu < min(*n,*il) || *iu > *n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
	if (*ldz < 1 || wantz && *ldz < *n) {
	    *info = -15;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = 1, i__2 = *n << 3;
	    if (*lwork < max(i__1,i__2) && ! lquery) {
		*info = -17;
	    }
	}
    }

    if (*info == 0) {
	nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__1, "DORMTR", uplo, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
	nb = max(i__1,i__2);
	lwkopt = (nb + 3) * *n;
	work[1] = (doublereal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYEVX", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	work[1] = 1.;
	return 0;
    }

    if (*n == 1) {
	work[1] = 7.;
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = a[a_dim1 + 1];
	} else {
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
		*m = 1;
		w[1] = a[a_dim1 + 1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.;
	}
	return 0;
    }

/*     Get machine constants. */

    safmin = dlamch_("Safe minimum", (ftnlen)12);
    eps = dlamch_("Precision", (ftnlen)9);
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
/* Computing MIN */
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
    rmax = min(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

    iscale = 0;
    abstll = *abstol;
    vll = *vl;
    vuu = *vu;
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	if (lower) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		dscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
/* L10: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		dscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
/* L20: */
	    }
	}
	if (*abstol > 0.) {
	    abstll = *abstol * sigma;
	}
	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

    indtau = 1;
    inde = indtau + *n;
    indd = inde + *n;
    indwrk = indd + *n;
    llwork = *lwork - indwrk + 1;
    dsytrd_(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwrk], &llwork, &iinfo, (ftnlen)1);
    lopt = (integer) (*n * 3 + work[indwrk]);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call DSTERF or DORGTR and SSTEQR.  If this fails for */
/*     some eigenvalue, then try DSTEBZ. */

    if ((alleig || indeig && *il == 1 && *iu == *n) && *abstol <= 0.) {
	dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
	indee = indwrk + (*n << 1);
	if (! wantz) {
	    i__1 = *n - 1;
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
	    dsterf_(n, &w[1], &work[indee], info);
	} else {
	    dlacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz, (
		    ftnlen)1);
	    dorgtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo, (ftnlen)1);
	    i__1 = *n - 1;
	    dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
	    dsteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[
		    indwrk], info, (ftnlen)1);
	    if (*info == 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ifail[i__] = 0;
/* L30: */
		}
	    }
	}
	if (*info == 0) {
	    *m = *n;
	    goto L40;
	}
	*info = 0;
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

    if (wantz) {
	*(unsigned char *)order = 'B';
    } else {
	*(unsigned char *)order = 'E';
    }
    indibl = 1;
    indisp = indibl + *n;
    indiwo = indisp + *n;
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwrk], &iwork[indiwo], info, (ftnlen)1, (ftnlen)1);

    if (wantz) {
	dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &
		ifail[1], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by DSTEIN. */

	indwkn = inde;
	llwrkn = *lwork - indwkn + 1;
	dormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

L40:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	dscal_(&imax, &d__1, &w[1], &c__1);
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

    if (wantz) {
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp1 = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp1) {
		    i__ = jj;
		    tmp1 = w[jj];
		}
/* L50: */
	    }

	    if (i__ != 0) {
		itmp1 = iwork[indibl + i__ - 1];
		w[i__] = w[j];
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
		if (*info != 0) {
		    itmp1 = ifail[i__];
		    ifail[i__] = ifail[j];
		    ifail[j] = itmp1;
		}
	    }
/* L60: */
	}
    }

/*     Set WORK(1) to optimal workspace size. */

    work[1] = (doublereal) lwkopt;

    return 0;

/*     End of DSYEVX */

} /* dsyevx_ */

