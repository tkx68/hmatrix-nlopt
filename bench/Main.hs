module Main where

import qualified AugLag as AugLag
import Criterion.Main
import qualified ISRES as ISRES

-- Our benchmark harness.
main :: IO ()
main = do
  let alkylX0 = [1.75, 1.2, 1.1, 3.048, 1.974, 0.893, 0.982, 8, 3.6, 1.45, 1, 1, 1, 1]
  defaultMain
    [ bgroup
        "Alkyl"
        [ bench "1" $ whnf ISRES.alkyl alkylX0,
          bench "1" $ whnf AugLag.alkyl alkylX0
        ],
      bgroup
        "ISRES01"
        [ bench "1" $ whnf ISRES.isres01 1,
          bench "10" $ whnf ISRES.isres01 10
        ],
      bgroup
        "ISRES02"
        [ bench "2" $ whnf ISRES.isres02 2,
          bench "10" $ whnf ISRES.isres02 10
        ],
      bgroup
        "ISRES03"
        [ bench "4" $ whnf ISRES.isres03 4,
          bench "10" $ whnf ISRES.isres03 10
        ]
    ]
