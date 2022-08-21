package rectfault

import "math"

func CalcDisp(a, x, y, z, c, dip, l, w, us, ud, ut float64) (ux, uy, uz float64, err error) {
	var ua, uh, ub, uc [3]float64

	// displacements
	// contribution of strike slip
	ua[0], ua[1], ua[2] = chinnerys(faStrike, a, x, y, z, c, dip, l, w)
	uh[0], uh[1], uh[2] = chinnerys(faStrike, a, x, y, -z, c, dip, l, w)
	ub[0], ub[1], ub[2] = chinnerys(fbStrike, a, x, y, z, c, dip, l, w)
	uc[0], uc[1], uc[2] = chinnerys(fcStrike, a, x, y, z, c, dip, l, w)
	uxs, uys, uzs := uxyz(ua, uh, ub, uc, z, dip)

	// contribution of dip slip
	ua[0], ua[1], ua[2] = chinnerys(faDip, a, x, y, z, c, dip, l, w)
	uh[0], uh[1], uh[2] = chinnerys(faDip, a, x, y, -z, c, dip, l, w)
	ub[0], ub[1], ub[2] = chinnerys(fbDip, a, x, y, z, c, dip, l, w)
	uc[0], uc[1], uc[2] = chinnerys(fcDip, a, x, y, z, c, dip, l, w)
	uxd, uyd, uzd := uxyz(ua, uh, ub, uc, z, dip)

	// contribution of tensile slip
	ua[0], ua[1], ua[2] = chinnerys(faTensile, a, x, y, z, c, dip, l, w)
	uh[0], uh[1], uh[2] = chinnerys(faTensile, a, x, y, -z, c, dip, l, w)
	ub[0], ub[1], ub[2] = chinnerys(fbTensile, a, x, y, z, c, dip, l, w)
	uc[0], uc[1], uc[2] = chinnerys(fcTensile, a, x, y, z, c, dip, l, w)
	uxt, uyt, uzt := uxyz(ua, uh, ub, uc, z, dip)

	ux = us*uxs + ud*uxd + ut*uxt
	uy = us*uys + ud*uyd + ut*uyt
	uz = us*uzs + ud*uzd + ut*uzt

	return ux, uy, uz, nil
}

// chinnerys performs Chinnery's operation with the given function f.
func chinnerys(f funcType, a, x, y, z, c, dip, l, w float64) (u1, u2, u3 float64) {
	cosd := math.Cos(dip)
	sind := math.Sin(dip)

	d := c - z
	p := y*cosd + d*sind
	q := y*sind - d*cosd

	var xis, etas [4]float64
	var v [4]okadaVars

	xis[0], etas[0] = x, p
	xis[1], etas[1] = x, p-w
	xis[2], etas[2] = x-l, p
	xis[3], etas[3] = x-l, p-w

	for i := 0; i < 4; i++ {
		xi, eta := xis[i], etas[i]
		v[i].R = math.Sqrt(xi*xi + eta*eta + q*q)
		v[i].Yt = eta*cosd + q*sind
		v[i].Dt = eta*sind - q*cosd
		v[i].Ct = v[i].Dt + z
		v[i].X = math.Sqrt(xi*xi + q*q)
		v[i].X11 = X11(xi, v[i].R)
		v[i].Y11 = Y11(eta, v[i].R)
		v[i].X32 = X32(xi, v[i].R)
		v[i].Y32 = Y32(eta, v[i].R)
		v[i].Z32 = Z32(eta, v[i].R, dip, q, z)

		// handling of singular values
		if math.Abs(q) < eps {
			v[i].Theta = 0.0
		} else {
			v[i].Theta = math.Atan((xi * eta) / (q * v[i].R))
		}

		// ln(R+xi)
		if math.Abs(v[i].R+xi) < eps {
			v[i].LnRpXi = -math.Log(v[i].R - xi)
		} else {
			v[i].LnRpXi = math.Log(v[i].R + xi)
		}

		// ln(R+eta)
		if math.Abs(v[i].R+eta) < eps {
			v[i].LnRpEta = -math.Log(v[i].R - eta)
		} else {
			v[i].LnRpEta = math.Log(v[i].R + eta)
		}
	}

	// perform Chinnery's operation:
	// f = f(x, p) - f(x, p-W), - f(x-L, p) + f(x-L, p-W)
	var f1, f2, f3 [4]float64
	for i := 0; i < 4; i++ {
		f1[i], f2[i], f3[i] = f(a, xis[i], etas[i], z, dip, q, v[i])
	}
	u1 = f1[0] - f1[1] - f1[2] + f1[3]
	u2 = f2[0] - f2[1] - f2[2] + f2[3]
	u3 = f3[0] - f3[1] - f3[2] + f3[3]

	return u1, u2, u3
}

// uxyz converts the disp (u1, u2, u3) to (ux, uy, uz).
func uxyz(ua, uh, ub, uc [3]float64, z, dip float64) (ux, uy, uz float64) {
	cosd := math.Cos(dip)
	sind := math.Sin(dip)

	ux = 1.0 / (2.0 * math.Pi) * (ua[0] - uh[0] + ub[0] + z*uc[0])
	uy = 1.0 / (2.0 * math.Pi) * ((ua[1]-uh[1]+ub[1]+z*uc[1])*cosd - (ua[2]-uh[2]+ub[2]+z*uc[2])*sind)
	uz = 1.0 / (2.0 * math.Pi) * ((ua[1]-uh[1]+ub[1]-z*uc[1])*sind + (ua[2]-uh[2]+ub[2]-z*uc[2])*cosd)

	return ux, uy, uz
}
