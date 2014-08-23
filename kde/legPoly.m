function val = legPoly(u, order)
% Returns the legendre polynomial of order 'order' at each element in u.
% u is an numPts x 1 vector
% This is sqrt( (2*order + 1)/ 2) *legendrePolynomial. The initial constant is
% needed for the kde

  switch order
    case 0
      val = ones(size(u));

    case 1
      val = u;

    case 2
      val = 1/2 * (3*u.^2 - 1);

    case 3
      val = 1/2 * (5*u.^3 - 3*u);

    case 4
      val = 1/8 * (35*u.^4 - 30*u.^2 + 3);

    otherwise
      temp = legendre(order, u);
      val = temp(1, :)';
      val = reshape(val, size(u));
  end

  % finally multiply by sqrt( (2*order + 1)/ 2)
  val = sqrt( (2*order + 1)/2 ) * val;
end

