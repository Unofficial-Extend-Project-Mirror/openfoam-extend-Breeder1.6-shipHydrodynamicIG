//scalarField period(2 * mathematicalConstant::pi / omega_);
// scalar factor(min(1.5 * period,db().time().value()) / (1.5 * period)); // Linear soft start
scalarField softStartFactor(Foam::sin(mathematicalConstant::pi / 2 * min(1.5 * wavePeriod_,db().time().value()) / (1.5 * wavePeriod_))); // Nonlinear soft start based on a sine

scalar eta(0.0);
// The current surface elevation at the boundary
eta = sum(ampl_ * Foam::cos(k_*distance - omega_ * db().time().value() + phase_) * softStartFactor) + seaLevel_;


