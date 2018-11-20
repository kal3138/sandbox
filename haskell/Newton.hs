module Newton (
    sqrt2,
    dsq2,
    newton
    ) where

sqrt2 :: Double -> Double
sqrt2 x = x^2 - 2

dsq2 :: Double -> Double
dsq2 x = 2*x

newton :: Double -> (Double -> Double) -> (Double -> Double) -> Double
newton x f df
    | abs((f x)/(df x)) < e = x
    | otherwise             = newton (x - (f x)/(df x)) f df
    where
        e = 1e-9 :: Double