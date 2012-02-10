harmonicNumbers = {};

function logfactorial (k)
{
	f = 0;
	for (; k>1; k = k-1)
	{
		f += Log(k);
	}
	return f;
}

function waringDensity (x, rho, phi)
{
	return Log(rho-1) + LnGamma (phi+rho) -LnGamma(phi+1) + LnGamma(x+phi) - LnGamma (x+phi+rho);
}

function yuleDensity (x, a)
{
	return Log (a) + Log (Beta(x, a+1));
}

function negBinomDensity (x, p, r)
{
	return LnGamma (r+x) - LnGamma(r) + r*Log(1-p) + x*Log(p) - logfactorial(x);
}

function paretoDensity (x, p)
{
	hsn = harmonicNumbers[p];
	if (hsn == 0)
	{
		vec = {1,10000} ["(_MATRIX_ELEMENT_COLUMN_+1)^(-p)"];
		hsn = +vec;
		harmonicNumbers[p] = +hsn;
	}
	return x^(-p)/hsn;
}

function likeFuncWaring (data, a, b)
{
	logL = 0;	
	data["sumupWaring"][""];
	return logL;
}

function likeFuncPareto (data, a)
{
	logL = 0;	
	data["sumupPareto"][""];
	return logL;
}

function likeFuncYule (data, rho)
{
	logL = 0;	
	data["sumupYule"][""];
	return logL;
}

function likeFuncNB (data, a, b)
{
	logL = 0;	
	data["sumupNB"][""];
	return logL;
}

function sumupWaring (key, value)
{
	logL += (0+value) * waringDensity (0 + key, a, b);
	return 0;
}

function sumupYule (key, value)
{
	logL += (0+value) * yuleDensity (0 + key, rho);
	return 0;
}

function sumupPareto (key, value)
{
	logL += (0+value) * Log(paretoDensity (0 + key, a));
	return 0;
}

function sumupNB (key, value)
{
	logL += value * negBinomDensity (0 + key, a, b);
	return 0;
}

totalNon0 = 0;	
maxD	  = 0;
degreeCounts = {};


for (k = 1; k < Abs(allDegs); k += 1)
{
	if (allDegs[k] > 0)
	{
		degreeCounts[k] += allDegs[k];
		totalNon0 += allDegs[k];
		maxD	   = Max (maxD, k);
	}
}


//fprintf (stdout, maxD, "\n", degreeCounts, "\n", totalNon0, "\n");
fprintf (stdout, "Degree,PDF,CDF\n");

sum = 0;
for (k = 1; k <= maxD; k+=1)
{
	if (degreeCounts[k]>0)
	{
		p = degreeCounts[k] / totalNon0;
		sum += p;
		
		fprintf (stdout, Join (",",{{k, p, sum}}), "\n");
	}
}

x2 :> -1;
x1 :> 1;
x1 :< 10000;
x2 :< 10000;

_degree_fit_results = {};

x1 = 3;
x2 = 0;

//VERBOSITY_LEVEL = 10;

Optimize (res, likeFuncWaring(degreeCounts,x1,x2));
fprintf (stdout, "\n\nWaring:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\nrho = ", x1, "\np = ", (x1-2)/(x1+x2-1), "\n");

_degree_fit_results ["Waring"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "rho": x1, "p": (x1-2)/(x1+x2-1)};

x2 := 0;
x1 :> 0;
Optimize (res, likeFuncYule(degreeCounts,x1));
fprintf (stdout, "\n\nYule:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\nrho = ", x1, "\n");

_degree_fit_results ["Yule"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "rho": x1};


x1 :< 1;
x1 :> 0;
x1 = 0.5;
x2 = 0.5;
x2 :> 0.000001;

Optimize (res, likeFuncNB(degreeCounts,x1,x2));
fprintf (stdout, "\n\nNeg. Binomial:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\np = ", x1, "\nr = ", x2, "\n");
_degree_fit_results ["Negative Binomial"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "p": x1};

x1 :> 0;
x1 :< 1e26;

/*for (z = 1; z < 3; z += 0.1)
{
	fprintf (stdout, z, ":", likeFuncPareto (degreeCounts, z), "\n");
}*/

Optimize (res, likeFuncPareto(degreeCounts,x1));
fprintf (stdout, "\n\nPareto:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\np = ", x1,"\n");
_degree_fit_results ["Pareto"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "rho": p};

function _THyPhyAskFor(key)
{
    if (key == "BestDistro")
    {
    	bestBIC = 1e100;
    	bestDistro = "";
        for (k = 0; k < Abs (_degree_fit_results); k += 1)
        {
        	key = _degree_fit_results["INDEXORDER"][k];
        	if ((_degree_fit_results[key])["BIC"] < bestBIC)
        	{
        		bestBIC = (_degree_fit_results[key])["BIC"];
        		bestDistro = key;
        	}
        }
        return bestDistro;
    }
    
    if (key == "Waring" || key == "Yule" || key == "Pareto")
    {
    	return (_degree_fit_results[key])["rho"];
    }
    
    if (key == "Negative Binomial")
    {
    	return 0.0;
    }
    return "_THyPhy_NOT_HANDLED_";
}
