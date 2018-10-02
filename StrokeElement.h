#ifndef STROKEELEMENT_H
#define STROKEELEMENT_H

#include <cmath>
#include <cfloat>
#include "colutils.h"

using namespace std;

struct StrokeElement
{
    int x, y; // position
    float u, v; // size (lengths of major and minor axes)
    float cos_theta, sin_theta; // cosine and sine of orientation
    float L, C, cos_H, sin_H; // color in lightness-chroma-hue

    const StrokeElement *nn[4]; // pointers to nearest neighbors in the four quadrants
    float nn_dist[4];

    StrokeElement()
    {
        do_nn_clear();
    }

    inline float get_orientation() const
    {
        return atan2(sin_theta, cos_theta);
    }
    inline void get_color(unsigned char & R, unsigned char & G, unsigned char & B) const
    {
        colutils::LabtoRGB(L, C * cos_H, C * sin_H, R, G, B);
    }

    void do_nn_clear()
    {
        for (int k = 0; k < 4; k++)
        {
            nn[k] = NULL;
            nn_dist[k] = FLT_MAX;
        }
    }
    void do_nn_check(const StrokeElement * const ptr_stroke)
    {
        static const float HALF_SQRT2 = sqrt(2) / 2;

        int dx = ptr_stroke->x - x;
        int dy = ptr_stroke->y - y;
        if (dx != 0 || dy != 0)
        {
            float radius = sqrt((float)(dx * dx + dy * dy));
            float cos_phi = (float)dx / radius;
            float sin_phi = (float)dy / radius;

            int quadrant;
            float cos_alpha = cos_phi * cos_theta + sin_phi * sin_theta;
            if (cos_alpha >= HALF_SQRT2)
                quadrant = 0;
            else if (cos_alpha <= - HALF_SQRT2)
                quadrant = 2;
            else if (sin_phi * cos_theta > cos_phi * sin_theta) // sin_alpha > 0
                quadrant = 1;
            else // sin_alpha < 0
                quadrant = 3;

            if (nn[quadrant] == NULL || nn_dist[quadrant] > radius)
            {
                nn[quadrant] = ptr_stroke;
                nn_dist[quadrant] = radius;
            }
        }
    }

    void rd_theta(const float lambda,
                  const float beta, const float prior_c, const float prior_s,
                  float &tmp_c, float &tmp_s) // one reaction-diffusion step
    {
        float sum_d = 0.0f, sum_w = 0.0f;
        for (int k = 0; k < 4; k++)
        {
            if (nn[k])
            {
                float w = 1.0f / nn_dist[k];
                sum_d += w * (nn[k]->sin_theta * cos_theta - nn[k]->cos_theta * sin_theta);
                sum_w += w;
            }
        }
        float d = lambda * (sum_d / sum_w) + beta * (prior_s * cos_theta - prior_c * sin_theta);
        float sin_d, cos_d;
        approx_sincos(d, sin_d, cos_d);
        tmp_c = cos_theta * cos_d - sin_theta * sin_d;
        tmp_s = sin_theta * cos_d + cos_theta * sin_d;
    }
    float rd_u(const float lambda,
               const float beta, const float prior,
               const float min_u, const float max_u) // one reaction-diffusion step
    {
        float sum_u = 0.0f, sum_w = 0.0f;
        for (int k = 0; k < 4; k++)
        {
            if (nn[k])
            {
                float w = 1.0f / nn_dist[k];
                sum_u += w * nn[k]->u;
                sum_w += w;
            }
        }
        return trim(u + lambda * (sum_u / sum_w - u) + beta * (prior - u), min_u, max_u);
    }
    float rd_v(const float lambda,
               const float beta, const float prior,
               const float min_v, const float max_v) // one reaction-diffusion step
    {
        float sum_v = 0.0f, sum_w = 0.0f;
        for (int k = 0; k < 4; k++)
        {
            if (nn[k])
            {
                float w = 1.0f / nn_dist[k];
                sum_v += w * nn[k]->v;
                sum_w += w;
            }
        }
        return trim(v + lambda * (sum_v / sum_w - v) + beta * (prior - v), min_v, max_v);
    }
    float rd_L(const float lambda,
               const float beta, const float prior,
               const float min_L = 0.0f, const float max_L = 100.0f) // one reaction-diffusion step
    {
        float sum_L = 0.0f, sum_w = 0.0f;
        for (int k = 0; k < 4; k++)
        {
            if (nn[k])
            {
                float w = 1.0f / nn_dist[k];
                sum_L += w * nn[k]->L;
                sum_w += w;
            }
        }
        return trim(L + lambda * (sum_L / sum_w - L) + beta * (prior - L), min_L, max_L);
    }
    float rd_C(const float lambda,
               const float beta, const float prior,
               const float min_C = 0.0f, const float max_C = 100.0f) // one reaction-diffusion step
    {
        float sum_C = 0.0f, sum_w = 0.0f;
        for (int k = 0; k < 4; k++)
        {
            if (nn[k])
            {
                float w = 1.0f / nn_dist[k];
                sum_C += w * nn[k]->C;
                sum_w += w;
            }
        }
        return trim(C + lambda * (sum_C / sum_w - C) + beta * (prior - C), min_C, max_C);
    }
    void rd_H(const float lambda,
              const float beta, const float prior_c, const float prior_s,
              float &tmp_c, float &tmp_s) // one reaction-diffusion step
    {
        float sum_d = 0.0f, sum_w = 0.0f;
        for (int k = 0; k < 4; k++)
        {
            if (nn[k])
            {
                float w = 1.0f / nn_dist[k];
                sum_d += w * (nn[k]->sin_H * cos_H - nn[k]->cos_H * sin_H);
                sum_w += w;
            }
        }
        float d = lambda * (sum_d / sum_w) + beta * (prior_s * cos_H - prior_c * sin_H);
        float sin_d, cos_d;
        approx_sincos(d, sin_d, cos_d);
        tmp_c = cos_H * cos_d - sin_H * sin_d;
        tmp_s = sin_H * cos_d + cos_H * sin_d;
    }

    void approx_sincos(const float rad, float & sin_rad, float & cos_rad)
    {
//        sin_x = sin(x);
//        cos_x = cos(x);

        static const float HALF_PI = 1.570796f;
        static const float C0 = sqrt(2) / 2;
        static const float C2 = 2.0f - 4.0f * C0;

        float x = rad / HALF_PI, tmp;

        while (x < -2) x += 4;
        while (x > 2) x -= 4;

        if (x < -1) // Quadrant III
        {
            x += 1.5f;
            tmp = C2 * x * x + C0;
            sin_rad = - tmp - x;
            cos_rad = - tmp + x;
        }
        else if (x < 0) // Quadrant IV
        {
            x += 0.5f;
            tmp = C2 * x * x + C0;
            sin_rad = - tmp + x;
            cos_rad =   tmp + x;
        }
        else if (x < 1) // Quadrant I
        {
            x -= 0.5f;
            tmp = C2 * x * x + C0;
            sin_rad =   tmp + x;
            cos_rad =   tmp - x;
        }
        else // Quadrant II
        {
            x -= 1.5f;
            tmp = C2 * x * x + C0;
            sin_rad =   tmp - x;
            cos_rad = - tmp - x;
        }
    }

    inline float trim(const float x, const float a, const float b)
    {
        assert(a < b);
        return x < a ? a : (x > b ? b : x);
    }
};

#endif // STROKEELEMENT_H
