package rectfault

import (
	"math"
)

// epsilon to check singularity
var eps = math.Nextafter(1.0, 2.0) - 1.0

// okadaVars stores common values required to compute displacements due to a rectangular source
type okadaVars struct {
	Xi, Eta                 float64
	R                       float64 // sqrt(xi^2+eta^2+p^2)
	Yt, Dt, Ct              float64 // y-tilde, d-tilde, c-tilde
	Theta                   float64 // atan(xi*eta/(q*R))
	X                       float64 // xi^2+q^2
	Cosd, Sind              float64 // cos(dip), sin(dip)
	LnRpXi                  float64 // dlog(R+xi)
	LnRpEta                 float64 // dlog(R+eta)
	X11, Y11, X32, Y32, Z32 float64
}

// funcType defines the types of functions to compute Chinney's operation.
// Following functions are available for Table6 of Okada (1992).
//
//	strike slip components: faStrike, fbStrike, fcStrike
//	dip slip components   : faDip, fbDip, fcDip
//	tensile components    : faTensile, fbTensile, fcTensile
type funcType func(float64, float64, float64, float64, float64, float64, okadaVars) (float64, float64, float64)

// functions that implement funcTypes
func faStrike(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	theta := v.Theta
	R := v.R
	lnRpEta := v.LnRpEta
	Y11 := v.Y11

	f1 = 0.5*theta + 0.5*a*xi*q*Y11
	f2 = 0.5 * a * q / R
	f3 = (1.0-a)/2.0*lnRpEta - (a/2.0)*(q*q)*Y11
	return f1, f2, f3
}

func fbStrike(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	theta := v.Theta
	R := v.R
	yt := v.Yt
	dt := v.Dt
	Y11 := v.Y11
	sind := v.Sind

	f1 = -xi*q*Y11 - theta - (1.0-a)/a*I1(xi, eta, dip, q, v)*sind
	f2 = -q/R + (1.0-a)/a*yt/(R+dt)*sind
	f3 = q*q*Y11 - (1.0-a)/a*I2(xi, eta, dip, q, v)*sind
	return f1, f2, f3
}

func fcStrike(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	R := v.R
	ct := v.Ct
	Y11 := v.Y11
	Z32 := v.Z32
	cosd := v.Cosd
	sind := v.Sind
	r3 := R * R * R

	f1 = (1.0-a)*xi*Y11*cosd - a*xi*q*Z32
	f2 = (1.0-a)*(cosd/R+2.0*q*Y11*sind) - a*ct*q/(r3)
	f3 = (1.0-a)*q*Y11*cosd - a*(ct*eta/(r3)-z*Y11+(xi*xi)*Z32)
	return f1, f2, f3
}

func faDip(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	theta := v.Theta
	R := v.R
	lnRpXi := v.LnRpXi
	X11 := v.X11

	f1 = a / 2.0 * q / R
	f2 = theta/2.0 + a/2.0*eta*q*X11
	f3 = (1-a)/2.0*lnRpXi - a/2.0*(q*q)*X11
	return f1, f2, f3
}

func fbDip(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	theta := v.Theta
	R := v.R
	dt := v.Dt
	X11 := v.X11
	cosd := v.Cosd
	sind := v.Sind

	f1 = -q/R + (1.0-a)/a*I3(xi, eta, dip, q, v)*sind*cosd
	f2 = -eta*q*X11 - theta - (1.0-a)/a*xi/(R+dt)*sind*cosd
	f3 = q*q*X11 + (1.0-a)/a*I4(xi, eta, dip, q, v)*sind*cosd
	return f1, f2, f3
}

func fcDip(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	R := v.R
	yt := v.Yt
	ct := v.Ct
	dt := v.Dt
	X11 := v.X11
	Y11 := v.Y11
	X32 := v.X32
	r3 := R * R * R
	cosd := v.Cosd
	sind := v.Sind

	f1 = (1.0-a)*cosd/R - q*Y11*sind - a*ct*q/r3
	f2 = (1.0-a)*yt*X11 - a*ct*eta*q*X32
	f3 = -dt*X11 - xi*Y11*sind - a*ct*(X11-q*q*X32)
	return f1, f2, f3
}

func faTensile(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	theta := v.Theta
	lnRpXi := v.LnRpXi
	lnRpEta := v.LnRpEta
	X11 := v.X11
	Y11 := v.Y11

	f1 = -(1.0-a)/2.0*lnRpEta - a/2.0*q*q*Y11
	f2 = -(1.0-a)/2.0*lnRpXi - a/2.0*q*q*X11
	f3 = theta/2.0 - a/2.0*q*(eta*X11+xi*Y11)
	return f1, f2, f3
}

func fbTensile(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	theta := v.Theta
	R := v.R
	dt := v.Dt
	X11 := v.X11
	Y11 := v.Y11
	sind := v.Sind

	f1 = q*q*Y11 - (1.0-a)/a*I3(xi, eta, dip, q, v)*(sind*sind)
	f2 = q*q*X11 + (1.0-a)/a*xi/(R+dt)*(sind*sind)
	f3 = q*(eta*X11+xi*Y11) - theta - (1.0-a)/a*I4(xi, eta, dip, q, v)*(sind*sind)
	return f1, f2, f3
}

func fcTensile(a, xi, eta, z, dip, q float64, v okadaVars) (f1, f2, f3 float64) {
	R := v.R
	dt := v.Dt
	yt := v.Yt
	ct := v.Ct
	X11 := v.X11
	Y11 := v.Y11
	X32 := v.X32
	Z32 := v.Z32
	cosd := v.Cosd
	sind := v.Sind

	f1 = -(1.0-a)*(sind/R+q*Y11*cosd) - a*(z*Y11-q*q*Z32)
	f2 = (1.0-a)*2.0*xi*Y11*sind + dt*X11 - a*ct*(X11-q*q*X32)
	f3 = (1.0-a)*(yt*X11+xi*Y11*cosd) + a*q*(ct*eta*X32+xi*Z32)
	return f1, f2, f3
}

// ---------- specific terms in Table 6 ----------
func I1(xi, eta, dip, q float64, v okadaVars) float64 {
	R := v.R
	dt := v.Dt

	return -xi/(R+dt)*math.Cos(dip) - I4(xi, eta, dip, q, v)*math.Sin(dip)
}

func I2(xi, eta, dip, q float64, v okadaVars) float64 {
	R := v.R
	dt := v.Dt
	return math.Log(R+dt) + I3(xi, eta, dip, q, v)*math.Sin(dip)
}

func I3(xi, eta, dip, q float64, v okadaVars) float64 {
	R := v.R
	dt := v.Dt
	yt := v.Yt
	lnRpEta := v.LnRpEta
	cosd := v.Cosd
	sind := v.Sind

	// avoid singularity
	if math.Abs(cosd) < eps {
		return 0.5 * (eta/(R+dt) + yt*q/((R+dt)*(R+dt)) - lnRpEta)
	}

	a := 1. / cosd
	return yt/(cosd*(R+dt)) - (a*a)*(lnRpEta-sind*math.Log(R+dt))
}

func I4(xi, eta, dip, q float64, v okadaVars) float64 {
	R := v.R
	X := v.X
	dt := v.Dt
	yt := v.Yt
	cosd := v.Cosd
	sind := v.Sind
	rt := R + dt

	// avoid singularity
	switch {
	case math.Abs(xi) < eps:
		return 0.
	case math.Abs(cosd) < eps:
		return 0.5 * (xi * yt) / (rt * rt)
	default:
		t1 := eta*(X+q*cosd) + X*(R+X)*sind
		t2 := xi * (R + X) * cosd

		return (sind / cosd * xi / rt) + (2.0/(cosd*cosd))*math.Atan(t1/t2)
	}
}

func X11(xi, R float64) float64 {
	if math.Abs(R+xi) < eps {
		return 0.
	}
	return 1. / (R * (R + xi))
}

func Y11(eta, R float64) float64 {
	if math.Abs(R+eta) < eps {
		return 0.
	}
	return 1.0 / (R * (R + eta))
}

func X32(xi, R float64) float64 {
	if math.Abs(R+xi) < eps {
		return 0.
	}

	r3 := R * R * R            // R^3
	rx2 := (R + xi) * (R + xi) // (R+xi)^2
	return (2.*R + xi) / (r3 * rx2)
}

func Y32(eta, R float64) float64 {
	if math.Abs(R+eta) < eps {
		return 0.
	}

	r3 := R * R * R              // R^3
	re2 := (R + eta) * (R + eta) // (R+xi)^2
	return (2.0*R + eta) / (r3 * re2)
}

func Z32(eta, R, dip, q, z float64) float64 {
	r3 := R * R * R // R^3
	sind := math.Sin(dip)
	return sind/r3 - h(dip, q, z)*Y32(eta, R)
}

func h(dip, q, z float64) float64 {
	cosd := math.Cos(dip)
	return q*cosd - z
}
