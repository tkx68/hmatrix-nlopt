module ISRES where

import Numeric.LinearAlgebra
import Numeric.NLOPT

isres01 :: Int -> Either Result Solution
isres01 n =
  let objf x = x `dot` x + 22 -- define objective
      stop = ObjectiveRelativeTolerance 1e-12 :| [] -- define stopping criterion
      algorithm = ISRES objf [] [] (SeedValue 22) -- specify algorithm
      lowerbounds = fromList $ replicate n (-10) -- specify bounds
      upperbounds = fromList $ replicate n (10) -- specify bounds
      problem = GlobalProblem lowerbounds upperbounds stop algorithm
      x0 = fromList $ replicate n 5 -- specify initial guess
   in minimizeGlobal problem x0

isres02 :: Int -> Either Result Solution
isres02 n =
  let objf x = - norm_1 x
      stop = ObjectiveRelativeTolerance 1e-12 :| [] -- define stopping criterion
      algorithm = ISRES objf [] [eqCon] (SeedValue 22) -- specify algorithm
      lowerbounds = fromList $ replicate n 0 -- specify bounds
      upperbounds = fromList $ replicate n 10 -- specify bounds
      eqCon =
        EqualityConstraint
          { eqConstraintFunctions = Scalar (\x -> x ! 1 * (x ! 0 - 5)),
            {-               Scalar
                            ( \x ->
                                let x1 = subVector 0 (n `div` 2) x
                                    x2 = subVector ((n `div` 2) + 1) n x
                                 in x1 `dot` (x2 - upperbounds)
                            ),
             -}
            eqConstraintTolerance = 1e-12
          }

      problem = GlobalProblem lowerbounds upperbounds stop algorithm
      x0 = fromList $ replicate n 5 -- specify initial guess
   in minimizeGlobal problem x0

isres03 :: Int -> Either Result Solution
isres03 n =
  if n < 3
    then error "n must be at least 3"
    else
      let objf x = - norm_1 x
          stop = ObjectiveRelativeTolerance 1e-12 :| [] -- define stopping criterion
          algorithm = ISRES objf [] [eqCon] (SeedValue 22) -- specify algorithm
          lowerbounds = fromList $ replicate n 0 -- specify bounds
          upperbounds = fromList $ replicate n 10 -- specify bounds
          eqCon =
            EqualityConstraint
              { eqConstraintFunctions =
                  Vector
                    2
                    ( \x _ ->
                        fromList
                          [ x ! 1 * (x ! 0 - 10),
                            x ! 2 * (x ! 1 - (10 - x ! 0))
                          ]
                    ),
                eqConstraintTolerance = 1e-12
              }

          problem = GlobalProblem lowerbounds upperbounds stop algorithm
          x0 = fromList $ replicate n 0 -- specify initial guess
       in minimizeGlobal problem x0

{-
See https://www.mat.univie.ac.at/~neum/glopt/coconut/Benchmark/Library1_new_v1.html

Alkyl problem in AMPL format:
var x2 := 1.745, >= 0, <= 2;
var x3 := 1.2, >= 0, <= 1.6;
var x4 := 1.1, >= 0, <= 1.2;
var x5 := 3.048, >= 0, <= 5;
var x6 := 1.974, >= 0, <= 2;
var x7 := 0.893, >= 0.85, <= 0.93;
var x8 := 0.928, >= 0.9, <= 0.95;
var x9 := 8, >= 3, <= 12;
var x10 := 3.6, >= 1.2, <= 4;
var x11 := 1.45, >= 1.45, <= 1.62;
var x12 := 1, >= 0.99, <= 1.01010101010101;
var x13 := 1, >= 0.99, <= 1.01010101010101;
var x14 := 1, >= 0.9, <= 1.11111111111111;
var x15 := 1, >= 0.99, <= 1.01010101010101;

minimize obj:  - 6.3*x5*x8 + 5.04*x2 + 0.35*x3 + x4 + 3.36*x6;

subject to
e2:  - 0.819672131147541*x2 + x5 - 0.819672131147541*x6 = 0;
e3: 0.98*x4 - x7*(0.01*x5*x10 + x4) = 0;
e4:  - x2*x9 + 10*x3 + x6 = 0;
e5: x5*x12 - x2*(1.12 + 0.13167*x9 - 0.0067*x9*x9) = 0;
e6: x8*x13 - 0.01*(1.098*x9 - 0.038*x9*x9) - 0.325*x7 - 0.57425 = 0;
e7: x10*x14 + 22.2*x11 - 35.82 = 0;
e8: x11*x15 - 3*x8 + 1.33 = 0;
-}
alkyl :: [Double] -> Either Result Solution
alkyl x0 =
  let lb = fromList [0, 0, 0, 0, 0, 0.85, 0.9, 3, 1.2, 1.45, 0.99, 0.99, 0.9, 0.99]
      ub = fromList [2, 1.6, 1.2, 5, 2, 0.93, 0.95, 12, 4, 1.62, 1.01010101010101, 1.01010101010101, 1.11111111111111, 1.01010101010101]
      objf x = - 6.3 * x ! 3 * x ! 6 + 5.04 * x ! 0 + 0.35 * x ! 1 + x ! 2 + 3.36 * x ! 4
      stop = ObjectiveRelativeTolerance 1e-12 :| []
      eqCons =
        [ EqualityConstraint
            { eqConstraintFunctions = Scalar (\x -> -0.819672131147541 * x ! 0 + x ! 3 - 0.819672131147541 * x ! 4),
              eqConstraintTolerance = 1e-10
            },
          EqualityConstraint
            { eqConstraintFunctions = Scalar (\x -> 0.98 * x ! 2 - x ! 5 * (0.01 * x ! 3 * x ! 8 + x ! 2)),
              eqConstraintTolerance = 1e-10
            },
          EqualityConstraint
            { eqConstraintFunctions = Scalar (\x -> - x ! 0 * x ! 7 + 10 * x ! 1 + x ! 4),
              eqConstraintTolerance = 1e-10
            },
          EqualityConstraint
            { eqConstraintFunctions = Scalar (\x -> x ! 3 * x ! 10 - x ! 0 * (1.12 + 0.13167 * x ! 7 - 0.0067 * x ! 7 * x ! 7)),
              eqConstraintTolerance = 1e-10
            },
          EqualityConstraint
            { eqConstraintFunctions = Scalar (\x -> x ! 6 * x ! 11 - 0.01 * (1.098 * x ! 7 - 0.038 * x ! 7 * x ! 7) - 0.325 * x ! 5 - 0.57425),
              eqConstraintTolerance = 1e-10
            },
          EqualityConstraint
            { eqConstraintFunctions = Scalar (\x -> x ! 8 * x ! 12 + 22.2 * x ! 9 - 35.82),
              eqConstraintTolerance = 1e-10
            },
          EqualityConstraint
            { eqConstraintFunctions = Scalar (\x -> x ! 9 * x ! 13 - 3 * x ! 6 + 1.33),
              eqConstraintTolerance = 1e-10
            }
        ]
      algorithm = ISRES objf [] eqCons (SeedValue 22)
      problem = GlobalProblem lb ub stop algorithm
   in minimizeGlobal problem (fromList x0)
