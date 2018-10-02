/**
  * Robust Gradient via Locally Linear Regression
  *
  * Copyright (C) 2010 Mingtian Zhao (mtzhao@stat.ucla.edu)
  *
  * Permission is hereby granted, free of charge, to any person
  * obtaining a copy of this software and associated documentation
  * files (the "Software"), to deal in the Software without
  * restriction, including without limitation the rights to use,
  * copy, modify, merge, publish, distribute, sublicense, and/or sell
  * copies of the Software, and to permit persons to whom the
  * Software is furnished to do so, subject to the following
  * conditions:
  *
  * The above copyright notice and this permission notice shall be
  * included in all copies or substantial portions of the Software.
  *
  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
  * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
  * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
  * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
  * OTHER DEALINGS IN THE SOFTWARE.
  */

#ifndef ROBUSTGRADIENT_VERSION
#define ROBUSTGRADIENT_VERSION 0

#define ROBUSTGRADIENT_M 4

#include <cassert>
#include <cmath>
#include <cstring>

using namespace std;

class RobustGradient
{
public:
    static void compute(const float image_data[], const bool mask_data[],
                        const int image_width, const int image_height,
                        float grad_x[], float grad_y[])
    {
        assert(image_data != NULL);
        assert(grad_x != NULL);
        assert(grad_y != NULL);

        float weights[ROBUSTGRADIENT_M + 1][ROBUSTGRADIENT_M + 1];
        for (int i = 0; i <= ROBUSTGRADIENT_M; i++)
        {
            for (int j = 0; j < i; j++)
            {
                weights[i][j] = weights[j][i] = exp(- (float)(i * i + j * j) / (float)(ROBUSTGRADIENT_M * ROBUSTGRADIENT_M));
            }
            weights[i][i] = exp(- (float)(i * i + i * i) / (float)(ROBUSTGRADIENT_M * ROBUSTGRADIENT_M));
        }

        for (int n = 0, y = 0; y < image_height; y++)
        {
            int min_yy = y > ROBUSTGRADIENT_M ? y - ROBUSTGRADIENT_M : 0;
            int max_yy = y + ROBUSTGRADIENT_M < image_height - 1 ? y + ROBUSTGRADIENT_M : image_height - 1;

            for (int x = 0; x < image_width; x++, n++)
            {
                if (mask_data[n])
                {
                    int min_xx = x > ROBUSTGRADIENT_M ? x - ROBUSTGRADIENT_M : 0;
                    int max_xx = x + ROBUSTGRADIENT_M < image_width - 1 ? x + ROBUSTGRADIENT_M : image_width - 1;

                    float a11 = 0.0f, a12 = 0.0f , a13 = 0.0f;
                    float a22 = 0.0f, a23 = 0.0f;
                    float a33 = 0.0f;
                    float b1 = 0.0f, b2 = 0.0f, b3 = 0.0f;

                    for (int yy = min_yy; yy <= max_yy; yy++)
                    {
                        const float *pz = image_data + min_xx + yy * image_width;
                        for (int xx = min_xx; xx <= max_xx; xx++)
                        {
                            int dx = xx - x;
                            int dy = yy - y;

                            float wt = weights[dx < 0 ? -dx : dx][dy < 0 ? -dy : dy];
                            float wx = wt * dx;
                            float wy = wt * dy;

                            a11 += wx * dx;
                            a12 += wy * dx;
                            a13 += wx;
                            a22 += wy * dy;
                            a23 += wy;
                            a33 += wt;

                            b1 += wx * (*pz);
                            b2 += wy * (*pz);
                            b3 += wt * (*pz);

                            pz++;
                        }
                    }

                    float det11 = a22 * a33 - a23 * a23;
                    float det12 = a12 * a33 - a13 * a23;
                    float det13 = a12 * a23 - a22 * a13;
                    float det22 = a11 * a33 - a13 * a13;
                    float det23 = a11 * a23 - a12 * a13;
                    float det = a11 * det11 - a12 * det12 + a13 * det13;

                    // Cramer's rule
                    grad_x[n] = (  b1 * det11 - b2 * det12 + b3 * det13) / det;
                    grad_y[n] = (- b1 * det12 + b2 * det22 - b3 * det23) / det;
                }
            }
        }
    }
};

#endif // ROBUSTGRADIENT_VERSION
