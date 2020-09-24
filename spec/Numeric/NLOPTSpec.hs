module Numeric.NLOPTSpec (spec) where

import Numeric.LinearAlgebra
import Numeric.NLOPT
import Test.Hspec

--
-- Next we solve the alkyl problem (see)
--
-- @
--     var x2 := 1.745, >= 0, <= 2;
--     var x3 := 1.2, >= 0, <= 1.6;
--     var x4 := 1.1, >= 0, <= 1.2;
--     var x5 := 3.048, >= 0, <= 5;
--     var x6 := 1.974, >= 0, <= 2;
--     var x7 := 0.893, >= 0.85, <= 0.93;
--     var x8 := 0.928, >= 0.9, <= 0.95;
--     var x9 := 8, >= 3, <= 12;
--     var x10 := 3.6, >= 1.2, <= 4;
--     var x11 := 1.45, >= 1.45, <= 1.62;
--     var x12 := 1, >= 0.99, <= 1.01010101010101;
--     var x13 := 1, >= 0.99, <= 1.01010101010101;
--     var x14 := 1, >= 0.9, <= 1.11111111111111;
--     var x15 := 1, >= 0.99, <= 1.01010101010101;

--     minimize obj:  - 6.3*x5*x8 + 5.04*x2 + 0.35*x3 + x4 + 3.36*x6;

--     subject to

--     e2:  - 0.819672131147541*x2 + x5 - 0.819672131147541*x6 = 0;
--     e3: 0.98*x4 - x7*(0.01*x5*x10 + x4) = 0;
--     e4:  - x2*x9 + 10*x3 + x6 = 0;
--     e5: x5*x12 - x2*(1.12 + 0.13167*x9 - 0.0067*x9*x9) = 0;
--     e6: x8*x13 - 0.01*(1.098*x9 - 0.038*x9*x9) - 0.325*x7 = 0.57425;
--     e7: x10*x14 + 22.2*x11 = 35.82;
--     e8: x11*x15 - 3*x8 = -1.33;
-- @
--
-- Expected solution:
--
-- @
--    x(1) = 3.03582
--    x(2) =  0.95
--    x(3) =  1.7037
--    x(4) =  0.543084
--    x(5) =  0.901319
--    x(6) =  10.4755
--    x(7) =  1.56164
--    x(8) =  1.53535
--    x(9) =  0.99
--    x(10) =  0.99
--    x(11) =  1.11111
--    x(12) =  0.99
--    x(13) =  1.58471
--    x(14) =  2
--    obj =  -1.765
-- @
spec :: Spec
spec = do
{-   describe "Sequential Least-Squares Quadratic Programming (SQSLP)" $ do
    it "solves the Alkyl problem" $ do
      let lb = fromList [0, 0, 0, 0, 0, 0.85, 0.9, 3, 1.2, 1.45, 0.99, 0.99, 0.9, 0.99]
          ub = fromList [2, 1.6, 1.2, 5, 2, 0.93, 0.95, 12, 4, 1.62, 1.01010101010101, 1.01010101010101, 1.11111111111111, 1.01010101010101]
          objfD x =
            ( -6.3 * x ! 3 * x ! 6 + 5.04 * x ! 0 + 0.35 * x ! 1 + x ! 2 + 3.36 * x ! 4,
              fromList [5.04, 0.35, 1, -6.3 * x ! 6, 3.36, 0, -6.3 * x ! 3, 0, 0, 0, 0, 0, 0, 0]
            )
          stop = ObjectiveRelativeTolerance 1e-8 :| []
          eqCons =
            [ EqualityConstraint
                { eqConstraintFunctions =
                    Vector
                      7
                      ( \x _ ->
                          (fromList
                            [ -0.819672131147541 * x ! 0 + x ! 3 - 0.819672131147541 * x ! 4,
                              0.98 * x ! 2 - x ! 5 * (0.01 * x ! 3 * x ! 8 + x ! 2),
                              - x ! 0 * x ! 7 + 10 * x ! 1 + x ! 4,
                              x ! 3 * x ! 10 - x ! 0 * (1.12 + 0.13167 * x ! 7 - 0.0067 * x ! 7 * x ! 7),
                              x ! 6 * x ! 11 - 0.01 * (1.098 * x ! 7 - 0.038 * x ! 7 * x ! 7) - 0.325 * x ! 5 - 0.57425,
                              x ! 8 * x ! 12 + 22.2 * x ! 9 - 35.82,
                              x ! 9 * x ! 13 - 3 * x ! 6 + 1.33
                            ],
                          toDense [
                            ((0,0), -0.819672131147541),
                            ((0,3), 1),
                            ((0,4), -0.819672131147541),
                            ((1,2), 0.98),
                            ((1,5), -(0.01 * x ! 3 * x ! 8 + x ! 2)),
                            ((1,3, -0.01*x ! 5 * x ! 8),
                            ((1,8, -0.01*x ! 5 * x ! 3),
                            ((1,2, -0.01*x ! 5),
                            ((2,0), -x ! 7),
                            ((2,0), -x ! 7),
                            ((2,7), -x ! 0),
                            ((2,1), 10),
                            ((2,4), 1),
                            ((3,3), x ! 10),
                            ((3,10), x ! 3),
                            ((3,0), -(1.12 + 0.13167 * x ! 7 - 0.0067 * x ! 7 * x ! 7)),
                            ((3,7), -x!0*(0.13167 -2*0.0067 * x ! 7)),
                            ((4,))
                            
                                ]
                      ),
                  eqConstraintTolerance = 1e-6
                }
            ]
          subproblem = LocalProblem (fromIntegral $ size lb) stop algorithm
          algorithm = SLSQP objfD [LowerBounds lb, UpperBounds ub] [{- InequalityConstraintsD -}] eqCons
          alkylX0 = [1.75, 1.2, 1.1, 3.048, 1.974, 0.893, 0.982, 8, 3.6, 1.45, 1, 1, 1, 1]
          res = minimizeLocal problem (fromList alkylX0)
      res
        `shouldBe` ( Right
                       ( Solution
                           { solutionCost = -1.765,
                             solutionParams =
                               fromList
                                 [ 3.03582,
                                   0.95,
                                   1.7037,
                                   0.543084,
                                   0.901319,
                                   10.4755,
                                   1.56164,
                                   1.53535,
                                   0.99,
                                   0.99,
                                   1.11111,
                                   0.99,
                                   1.58471,
                                   2
                                 ],
                             solutionResult = FTOL_REACHED
                           }
                       )
                   )
 -}
  describe "Augmented Lagrangian based on SBPLX (Subplex)" $ do
    it "solves the Alkyl problem" $ do
      let lb = fromList [0, 0, 0, 0, 0, 0.85, 0.9, 3, 1.2, 1.45, 0.99, 0.99, 0.9, 0.99]
          ub = fromList [2, 1.6, 1.2, 5, 2, 0.93, 0.95, 12, 4, 1.62, 1.01010101010101, 1.01010101010101, 1.11111111111111, 1.01010101010101]
          objf x = -6.3 * x ! 3 * x ! 6 + 5.04 * x ! 0 + 0.35 * x ! 1 + x ! 2 + 3.36 * x ! 4
          stop = ObjectiveRelativeTolerance 1e-8 :| []
          eqCons =
            [ EqualityConstraint
                { eqConstraintFunctions =
                    Vector
                      7
                      ( \x _ ->
                          fromList
                            [ -0.819672131147541 * x ! 0 + x ! 3 - 0.819672131147541 * x ! 4,
                              0.98 * x ! 2 - x ! 5 * (0.01 * x ! 3 * x ! 8 + x ! 2),
                              - x ! 0 * x ! 7 + 10 * x ! 1 + x ! 4,
                              x ! 3 * x ! 10 - x ! 0 * (1.12 + 0.13167 * x ! 7 - 0.0067 * x ! 7 * x ! 7),
                              x ! 6 * x ! 11 - 0.01 * (1.098 * x ! 7 - 0.038 * x ! 7 * x ! 7) - 0.325 * x ! 5 - 0.57425,
                              x ! 8 * x ! 12 + 22.2 * x ! 9 - 35.82,
                              x ! 9 * x ! 13 - 3 * x ! 6 + 1.33
                            ]
                      ),
                  eqConstraintTolerance = 1e-6
                }
            ]
          problem = AugLagProblem eqCons [] (AUGLAG_EQ_LOCAL subproblem)
          subproblem = LocalProblem (fromIntegral $ size lb) stop algorithm
          algorithm = SBPLX objf [LowerBounds lb, UpperBounds ub] Nothing
          alkylX0 = [1.75, 1.2, 1.1, 3.048, 1.974, 0.893, 0.982, 8, 3.6, 1.45, 1, 1, 1, 1]
          res = minimizeAugLag problem (fromList alkylX0)
      res
        `shouldBe` ( Right
                       ( Solution
                           { solutionCost = -1.765,
                             solutionParams =
                               fromList
                                 [ 3.03582,
                                   0.95,
                                   1.7037,
                                   0.543084,
                                   0.901319,
                                   10.4755,
                                   1.56164,
                                   1.53535,
                                   0.99,
                                   0.99,
                                   1.11111,
                                   0.99,
                                   1.58471,
                                   2
                                 ],
                             solutionResult = FTOL_REACHED
                           }
                       )
                   )
