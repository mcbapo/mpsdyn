double ramp_xvalue(double x_start,double x_end,int steps_total,int step)
{
  double result;
  //First option: linear ramping
  //result = x_start + (x_end-x_start)*((double) step)/((double) steps_total);
  
  //Second option: exponentially fast to the end
  //result = x_end*(exp((double) step)-1.0)/(exp((double) steps_total)-1.0);
  
  //Third option as polynomial, the order of the polynomial decides how flat the curve is
  result = (x_end-x_start)*pow((double) step,3.0)/pow((double) steps_total,3.0) + x_start;
  //This version is needed if I want to read a ground state which is already at a finite value of x
  //result = x_end*pow((double) step,3.0)/pow((double) steps_total,3.0);
  
  return result;
}

double ramp_xvalue_mixed(double x_start,double x_end,int steps_total,int step,int tchange)
{
  double result,a,b,c,d;
  
  a = -1.0*(pow(((double)steps_total),3.0)*(((double)steps_total)*x_start-((double)tchange)*x_start-x_end+x_start))/((3.0*((double)steps_total)-2.0*((double)tchange))*pow(((double)tchange),2.0));
  b = x_start;
  c = (((double)tchange)*x_start+3.0*x_end-3.0*x_start)/(3.0*((double)steps_total)-2.0*((double)tchange));
  d = -1.0*(((double)tchange)*((double)steps_total)*x_start-3.0*((double)steps_total)*x_start+2.0*((double)tchange)*x_end)/(3.0*((double)steps_total)-2.0*((double)tchange));
  
  if(step <= tchange)
    result= a*pow(((double)step)/((double)steps_total),3.0)+b;
  else
    result = c*((double)step)+d;
  
  return result;
}


double ramp_xvalue(double x_start,double x_end,int steps_total,int step, int offset)
{
  double result;
  //First option: linear ramping
  //result = x_start + (x_end-x_start)*((double) step)/((double) steps_total);
  
  //Second option: exponentially fast to the end
  //result = x_end*(exp((double) step)-1.0)/(exp((double) steps_total)-1.0);
  
  //Third option as polynomial, the order of the polynomial decides how flat the curve is
  result = x_end*pow((double) (step+offset),3.0)/pow((double) steps_total,3.0);
  
  return result;
}


//Ramp x according to a polynomial in time, the order is given by order
double ramp_xvalue_time(double x_start, double x_end, double total_time, double elapsed_time,double order)
{
  double result;
  result = (x_end-x_start)*pow(elapsed_time,order)/pow(total_time,order) + x_start;
  
  return result;
}

//Ramp x according to a third order polynomial in time
//Function is more ore less deprecated, I just keep it as wrapper to the new version, to have compatibility to my old codes
double ramp_xvalue_time(double x_start, double x_end, double total_time, double elapsed_time)
{
  return ramp_xvalue_time(x_start,x_end,total_time,elapsed_time,3.0);
}

//Idea, first I follow the simple curve x(t)=x_end*(t/T_1)^3. At some time ts where I reached a value xs, I switch to another, faster curve such that everything is done in T2<T1. The subsequent curve is smoothly connected to the first one
//Attention, the switching time has to be designed such that ts<T2, otherwise the thing goes crazy
//Still experimental!!!
double ramp_xvalue_time_hybrid(double x_start, double x_end, double total_time, double elapsed_time, double T1)
{
  double result;
  
  //Up to now it works only for x_start=0, therefore check for that
  if(x_start != 0)
	cout << "Warning, x_start should be 0, current value is " << x_start << endl;
  
  //Compute switching time
  double ts = pow(10.0/x_end,1.0/3.0)*T1;
  
  if(ts>total_time)
  {
	  cout << "Error, switching time is larger than total time" << endl;
	  cout << "Program will be aborted" << endl;
	  exit(666);
  }
  
  if(elapsed_time <= ts)
	result = x_end*pow(elapsed_time,3.0)/pow(T1,3.0);
  else
  {
	  double nenner=pow(T1,3.0)*(pow(total_time,3.0)-3.0*total_time*ts*ts+2.0*pow(ts,3.0));
	  double b=x_end*(pow(T1,3.0)-3.0*total_time*ts*ts+2.0*pow(ts,3.0))/nenner;
	  double d=-3.0*x_end*ts*ts*(pow(T1,3.0)-pow(total_time,3.0))/nenner;
	  double e=2.0*x_end*pow(ts,3.0)*(pow(T1,3.0)-pow(total_time,3.0))/nenner;
	  
	  result = b*pow(elapsed_time,3.0)+d*elapsed_time+e;
  }
  
  return result;
}

double ramp_muvalue(double mu_start,double mu_end,int steps_total, int step)
{
  double result;
  //First option: linear decrease/increase
  //result = mu_start + (mu_end-mu_start)*((double) step)/((double) steps_total);
  
  //Option as third order polynomial
  result = (mu_end-mu_start)*pow((double) step,3.0)/pow((double) steps_total,3.0) + mu_start;
   
  return result;
}

//Basically the same as the first function in the file, but this time I only have the linear option, such that I can directly choose which ramping is chosen by calling the right function and I don't have to comment / uncomment the right option all the time
double ramp_xvalue_linear(double x_start,double x_end,int steps_total,int step)
{
  double result;
  //linear ramping
  result = x_start + (x_end-x_start)*((double) step)/((double) steps_total);
  return result;
}

double ramp_xvalue_cubic(double x_start,double x_end,int steps_total,int step)
{
  double result;  
  //Third option as third order polynomial
  result = (x_end-x_start)*pow((double) step,3.0)/pow((double) steps_total,3.0) + x_start;
  return result;
}


//Function to ramp the BOND DIMENSION during evolution (on a long run I should move this one to a more sensible place, however for a quick test, I place it here)
// int ramp_bond(int D_start,int D_end,int steps_total,int step)
// {
//   return (int)((D_end-D_start)*(((double)step)/((double) steps_total))+D_start);
// }

int ramp_bond(int D_start, int D_end, int steps_total, int step)
{
  int result;
  result = (int)((D_end-D_start)*pow((double)step,3.0)/pow((double)steps_total,3.0)) + D_start;
  
  return result;
}
