%% Main function to generate tests
function tests = interpolationTest
  tests = functiontests(localfunctions);
end

%% Test Functions
function interp1cubicLTest(testCase)
  % Test specific code
  % This should be linear in [1,2] and [3,4] and equal to f in [2,3]
  tol = 1e-10;
  % Knots
  x = 1:4;
  f = @(x) 3*x.^3 + 4*x.^2 + x + 1;
  % Linear in [x(1),x(2)]
  y1 = x(1):0.01:x(2);
  out1 = interp1cubicL(x,f(x),y1);
  expected1 = f(x(1)) + (f(x(2))-f(x(1)))/(x(2)-x(1))*(y1 - x(1));
  res1 = norm(out1 - expected1);
  
  assert(res1<tol)
  % Equal to f in [x(2),x(3)]
  y2 = x(2):0.01:x(3);
  out2 = interp1cubicL(x,f(x),y2);
  expected2 = f(y2);
  res2 = norm(out2 - expected2);
  
  assert(res2<tol)
  % Linear in [x(3),x(4)]
  y3 = x(3):0.01:x(4);
  out3 = interp1cubicL(x,f(x),y3);
  expected3 = f(x(3)) + (f(x(4))-f(x(3)))/(x(4)-x(3))*(y3 - x(3));
  res3 = norm(out3 - expected3);
  
  assert(res3<tol)
end % function interp1cubicLTest

function interp3limTest(testCase)
  % Test to see that interp3lim in locally monotonic
  X = (0:6);
  M = sin(X);
  y = 0:0.1:6;
  M_y = interp3lim(X,M,y);
  assert(all(and(M_y<=max(M),M_y>=min(M))));
end % function interp3limTest

%
%
%function interp_ENO_coefsTest(testCase)
%end
%
%function interp_ENOTest(testCase)
%end
%
%function ppval_limTest(testCase)
%end
