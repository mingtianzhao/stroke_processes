/**
  * Color Utilitiy Toolkit
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

#ifndef COLUTILS_VERSION
#define COLUTILS_VERSION 0

struct colutils
{
    typedef unsigned char uchar;

    inline static void RGBtoLab(const uchar R, const uchar G, const uchar B, float &L, float &a, float &b)
    {
        float X, Y, Z;
        colutils::RGBtoXYZ(R, G, B, X, Y, Z);
        colutils::XYZtoLab(X, Y, Z, L, a, b);
    }
    inline static void LabtoRGB(const float L, const float a, const float b, uchar &R, uchar &G, uchar &B)
    {
        float X, Y, Z;
        colutils::LabtoXYZ(L, a, b, X, Y, Z);
        colutils::XYZtoRGB(X, Y, Z, R, G, B);
    }

    static void RGBtoXYZ(const uchar R, const uchar G, const uchar B, float &X, float &Y, float &Z)
    {
        static const float M[3][3] =
        {
            { 0.4124f / 255.0f, 0.3576f / 255.0f, 0.1805f / 255.0f },
            { 0.2126f / 255.0f, 0.7152f / 255.0f, 0.0722f / 255.0f },
            { 0.0193f / 255.0f, 0.1192f / 255.0f, 0.9505f / 255.0f }
        };

        X = M[0][0] * (float)R + M[0][1] * (float)G + M[0][2] * (float)B;
        Y = M[1][0] * (float)R + M[1][1] * (float)G + M[1][2] * (float)B;
        Z = M[2][0] * (float)R + M[2][1] * (float)G + M[2][2] * (float)B;
    }
    static void XYZtoRGB(const float X, const float Y, const float Z, uchar &R, uchar &G, uchar &B)
    {
        static const float M[3][3] =
        {
            { 3.2410f * 255.0f,-1.5374f * 255.0f,-0.4986f * 255.0f },
            {-0.9692f * 255.0f, 1.8760f * 255.0f, 0.0416f * 255.0f },
            { 0.0556f * 255.0f,-0.2040f * 255.0f, 1.0570f * 255.0f }
        };

        float fR = M[0][0] * X + M[0][1] * Y + M[0][2] * Z;
        float fG = M[1][0] * X + M[1][1] * Y + M[1][2] * Z;
        float fB = M[2][0] * X + M[2][1] * Y + M[2][2] * Z;

        R = (fR < 0) ? 0 : ((fR > 255) ? 255 : (uchar)fR);
        G = (fG < 0) ? 0 : ((fG > 255) ? 255 : (uchar)fG);
        B = (fB < 0) ? 0 : ((fB > 255) ? 255 : (uchar)fB);
    }

    static void XYZtoLab(float X, float Y, float Z, float &L, float &a, float &b)
    {
        static const float D65[3] = { 0.9505f, 1.0000f, 1.0890f };
        static const float delta = (6.0f / 29.0f);
        static const float t0 = delta * delta * delta;
        static const float C1 = 1.0f / (3.0f * delta * delta);
        static const float C0 = delta * (2.0f / 3.0f);

        X /= D65[0];
        Y /= D65[1];
        Z /= D65[2];

        X = (X > t0) ? approx_cbrt(X) : (C0 + C1 * X);
        Y = (Y > t0) ? approx_cbrt(Y) : (C0 + C1 * Y);
        Z = (Z > t0) ? approx_cbrt(Z) : (C0 + C1 * Z);

        L = 116.0f * Y - 16.0f;
        a = 500.0f * (X - Y);
        b = 200.0f * (Y - Z);
    }
    static void LabtoXYZ(const float L, const float a, const float b, float &X, float &Y, float &Z)
    {
        static const float D65[3] = { 0.9505f, 1.0000f, 1.0890f };
        static const float delta = (6.0f / 29.0f);
        static const float t = (3.0f * delta * delta);
        static const float f0 = (16.0f / 116.0f);

        float fy = (L + 16.0f) / 116.0f;
        float fx = fy + a / 500.0f;
        float fz = fy - b / 200.0f;

        X = (fx > delta) ? (D65[0] * fx * fx * fx) : ((fx - f0) * t * D65[0]);
        Y = (fy > delta) ? (D65[1] * fy * fy * fy) : ((fy - f0) * t * D65[1]);
        Z = (fz > delta) ? (D65[2] * fz * fz * fz) : ((fz - f0) * t * D65[2]);
    }

    static float approx_cbrt(float x)
    {
        if (x < 0) return -approx_cbrt(-x);

        float y = x;
        unsigned char *p = (unsigned char *)(&y);
        *p = *p / 3 + 709921077;

        float t = y * y * y;
        y *= (t + x + x) / (t + t + x);
        t = y * y * y;
        y *= (t + x + x) / (t + t + x);

        return y;
    }
};

#endif // COLUTILS_VERSION
