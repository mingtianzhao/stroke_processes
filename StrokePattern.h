#ifndef STROKEPATTERN_H
#define STROKEPATTERN_H

#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include "CImg.h"
#include "RobustGradient.h"
#include "colutils.h"
#include "StrokeElement.h"

using namespace std;
using namespace cimg_library;

class StrokePattern
{
private:
    struct SAL_ordering
    {
        bool operator ()(const pair<int, float> & a, const pair<int, float> & b) const
        {
            return a.second > b.second;
        }
    };

public:
    StrokePattern(const CImg<unsigned char> & src,
                  const CImg<float> & orientation, // orientation map in radians
                  const CImg<bool> & msk)
    {
        width = src.width();
        height = src.height();
        num_pixels = width * height;
        num_valid_pixels = 0;
        for (int n = 0; n < num_pixels; n++)
        {
            if (msk[n])
                num_valid_pixels++;
        }

        reseed();

        const unsigned char * const pR = src.data();
        const unsigned char * const pG = pR + num_pixels;
        const unsigned char * const pB = pG + num_pixels;

        refLIT.assign(width, height, 1, 1, 0.0f);
        refCHR.assign(width, height, 1, 1, 0.0f);
        refC_H.assign(width, height, 1, 1, 0.0f);
        refS_H.assign(width, height, 1, 1, 0.0f);

        for (int n = 0; n < num_pixels; n++)
        {
            if (msk[n])
            {
                colutils::RGBtoLab(pR[n], pG[n], pB[n], refLIT[n], refC_H[n], refS_H[n]);
                refCHR[n] = sqrt(refC_H[n] * refC_H[n] + refS_H[n] * refS_H[n]);
                if (refCHR[n] > 1e-9f)
                {
                    refC_H[n] /= refCHR[n];
                    refS_H[n] /= refCHR[n];
                }
                else
                {
                    refC_H[n] = 1.0f;
                    refS_H[n] = 0.0f;
                }
            }
        }

        refG_X.assign(width, height, 1, 1, 0.0f);
        refG_Y.assign(width, height, 1, 1, 0.0f);

        RobustGradient::compute(refLIT.data(), msk.data(), width, height, refG_X.data(), refG_Y.data());

        refSAL.assign(width, height, 1, 1, 0.0f);

        float max_sal = 0.0f;
        vector< pair<int, float> > vecSAL;
        vecSAL.reserve(num_valid_pixels);
        for (int n = 0; n < num_pixels; n++)
        {
            if (msk[n])
            {
                refSAL[n] = sqrt(refG_X[n] * refG_X[n] + refG_Y[n] * refG_Y[n]);
                refG_X[n] =-sin(orientation[n]);
                refG_Y[n] = cos(orientation[n]);

                if (max_sal < refSAL[n])
                    max_sal = refSAL[n];

                vecSAL.push_back(make_pair(n, refSAL[n]));
            }
        }
        refSAL /= max_sal;
        sort(vecSAL.begin(), vecSAL.end(), SAL_ordering());

        refSIZ.assign(width, height, 1, 1, 0.0f);
        static const float rateSAL = 2.5f;
        static const float rateSIZ = 7.5f;
        for (int n = 0; n < num_valid_pixels; n++)
        {
            float t = (float)n / (float)(num_valid_pixels - 1);
            refSAL[vecSAL[n].first] = - log(1.0f / exp(rateSAL) + t * (1.0f - 1.0f / exp(rateSAL))) / rateSAL;
            if (t > 0.5f)
                refSIZ[vecSAL[n].first] = 0.5f - 0.5f * log(1.0f / exp(rateSIZ) + 2.0f * (1.0f - t) * (1.0f - 1.0f / exp(rateSIZ))) / rateSIZ;
            else if (t < 0.5f)
                refSIZ[vecSAL[n].first] = 0.5f + 0.5f * log(1.0f / exp(rateSIZ) + 2.0f * (       t) * (1.0f - 1.0f / exp(rateSIZ))) / rateSIZ;
        }
    }

    StrokePattern(const CImg<unsigned char> & src,
                  const CImg<bool> & msk)
    {
        width = src.width();
        height = src.height();
        num_pixels = width * height;
        num_valid_pixels = 0;
        for (int n = 0; n < num_pixels; n++)
        {
            if (msk[n])
                num_valid_pixels++;
        }

        reseed();

        const unsigned char * const pR = src.data();
        const unsigned char * const pG = pR + num_pixels;
        const unsigned char * const pB = pG + num_pixels;

        refLIT.assign(width, height, 1, 1, 0.0f);
        refCHR.assign(width, height, 1, 1, 0.0f);
        refC_H.assign(width, height, 1, 1, 0.0f);
        refS_H.assign(width, height, 1, 1, 0.0f);

        for (int n = 0; n < num_pixels; n++)
        {
            if (msk[n])
            {
                colutils::RGBtoLab(pR[n], pG[n], pB[n], refLIT[n], refC_H[n], refS_H[n]);
                refCHR[n] = sqrt(refC_H[n] * refC_H[n] + refS_H[n] * refS_H[n]);
                if (refCHR[n] > 1e-9f)
                {
                    refC_H[n] /= refCHR[n];
                    refS_H[n] /= refCHR[n];
                }
                else
                {
                    refC_H[n] = 1.0f;
                    refS_H[n] = 0.0f;
                }
            }
        }

        refG_X.assign(width, height, 1, 1, 0.0f);
        refG_Y.assign(width, height, 1, 1, 0.0f);

        RobustGradient::compute(refLIT.data(), msk.data(), width, height, refG_X.data(), refG_Y.data());

        refSAL.assign(width, height, 1, 1, 0.0f);

        float max_sal = 0.0f;
        vector< pair<int, float> > vecSAL;
        vecSAL.reserve(num_valid_pixels);
        for (int n = 0; n < num_pixels; n++)
        {
            if (msk[n])
            {
                refSAL[n] = sqrt(refG_X[n] * refG_X[n] + refG_Y[n] * refG_Y[n]);
                if (refSAL[n] > 1e-9f)
                {
                    refG_X[n] /= refSAL[n];
                    refG_Y[n] /= refSAL[n];
                }
                else
                {
                    refG_X[n] = 1.0f;
                    refG_Y[n] = 0.0f;
                }

                if (max_sal < refSAL[n])
                    max_sal = refSAL[n];

                vecSAL.push_back(make_pair(n, refSAL[n]));
            }
        }
        refSAL /= max_sal;
        sort(vecSAL.begin(), vecSAL.end(), SAL_ordering());

        refSIZ.assign(width, height, 1, 1, 0.0f);
        static const float rateSAL = 2.5f;
        static const float rateSIZ = 7.5f;
        for (int n = 0; n < num_valid_pixels; n++)
        {
            float t = (float)n / (float)(num_valid_pixels - 1);
            refSAL[vecSAL[n].first] = - log(1.0f / exp(rateSAL) + t * (1.0f - 1.0f / exp(rateSAL))) / rateSAL;
            if (t > 0.5f)
                refSIZ[vecSAL[n].first] = 0.5f - 0.5f * log(1.0f / exp(rateSIZ) + 2.0f * (1.0f - t) * (1.0f - 1.0f / exp(rateSIZ))) / rateSIZ;
            else if (t < 0.5f)
                refSIZ[vecSAL[n].first] = 0.5f + 0.5f * log(1.0f / exp(rateSIZ) + 2.0f * (       t) * (1.0f - 1.0f / exp(rateSIZ))) / rateSIZ;
        }
    }

    void reseed()
    {
        rand_seed = (unsigned int)time(NULL);

        num_strokes = 0;

        ctvDEN = -1;
        ctvORD = -1;
        ctvLEN = -1;
        ctvTHK = -1;
        ctvLIT = -1;
        ctvCHR = -1;
        ctvHUE = -1;
    }

    void sample(const int density,
                const int non_uniformity,
                const int local_isotropy,
                const int coarseness,
                const int size_contrast,
                const int lightness_contrast,
                const int chroma_contrast,
                const int hue_contrast)
    {
        reaction_diffusion_cycles = 5;

        compute_pos(5 + density * 5, non_uniformity, 0.1f + coarseness * 0.1f);
        compute_ori(100 - local_isotropy);

        compute_size_u(1.0f + coarseness * 0.5f, 2.5f + coarseness * 1.5f, size_contrast);
        compute_size_v(0.2f + coarseness * 0.2f, 0.6f + coarseness * 0.6f, size_contrast);

        compute_color_L(lightness_contrast);
        compute_color_C(chroma_contrast);
        compute_color_H(hue_contrast);
    }

    inline int get_width() const { return width; }
    inline int get_height() const { return height; }
    inline int get_num_strokes() const { return num_strokes; }
    inline const vector<StrokeElement> & get_strokes() const { return strokes; }
    inline const StrokeElement & get_stroke(const int i) const { return strokes[i]; }

private:
    unsigned int rand_seed;
    int width, height, num_pixels, num_valid_pixels, num_strokes;
    int reaction_diffusion_cycles;

    float minLEN, maxLEN, minTHK, maxTHK; // stroke size limits

    // contrast values and corresponding reference maps
    int ctvDEN, ctvORD; // density and orientation
    int ctvLEN, ctvTHK; // length and thickness
    int ctvLIT, ctvCHR, ctvHUE; // lightness, chroma and hue
    CImg<float> refLIT, refCHR, refC_H, refS_H, refG_X, refG_Y, refSAL, refSIZ;

    CImg<float> pos_cmf;

    vector<StrokeElement> strokes;

private:
    void compute_pos(const int stroke_count, const int contrast, const float inhibition_radius)
    {
        bool resample = false;

        if (num_strokes != stroke_count)
        {
            num_strokes = stroke_count;
            resample = true;
        }

        if (ctvDEN != contrast)
        {
            ctvDEN = contrast;
            resample = true;

            const float power = 1.0f + (float)(contrast - 50) / 100.0f;

            pos_cmf = refSAL;
            pos_cmf[0] = pow(pos_cmf[0], power);
            for (int n = 1; n < num_pixels; n++)
            {
                pos_cmf[n] = pos_cmf[n - 1] + pow(pos_cmf[n], power);
            }
            pos_cmf /= pos_cmf[num_pixels - 1];
        }

        if (resample)
        {
            CImg<bool> inhibition_mask(width, height, 1, 1, true);
            bool b_false = false;

            srand(rand_seed);
            strokes.clear();
            strokes.reserve(num_strokes);
            for (int i = 0; i < num_strokes; i++)
            {
                int a, b, num_trials = 0;
                do {
                    float tmp = (float)rand() / (float)(RAND_MAX + 1); // in [0,1)

                    // binary search
                    a = 0;
                    b = num_pixels - 1;
                    while (a != b)
                    {
                        int n = (a + b) >> 1;
                        if (tmp < pos_cmf[n]) b = n;
                        else a = n + 1;
                    };

                    num_trials++;
                } while (!inhibition_mask[a] && num_trials < 20);

                if (inhibition_mask[a])
                {
                    div_t divresult = div(a, width);
                    StrokeElement se;
                    se.x = divresult.rem;
                    se.y = divresult.quot;
                    strokes.push_back(se);

                    inhibition_mask.draw_circle(se.x, se.y, inhibition_radius, &b_false);
                }
            }
            num_strokes = (int)strokes.size();

            // must redo compute_ori
            ctvORD = -1;
        }
    }

    void compute_ori(const int contrast)
    {
        if (ctvORD != contrast)
        {
            ctvORD = contrast;

            // assign orientation values according to the orientation field
            // with PI/2 rotation from gradient orientation
            for (int i = 0; i < num_strokes; i++)
            {
                const int n = stroke_offset(i);
                strokes[i].cos_theta = refG_Y[n];
                strokes[i].sin_theta =-refG_X[n];
            }

            // compute graph topology according to the orientation field
            compute_graph_topology();

            if (contrast != 50)
            {
                const float neighbor_weight = (float)(50 - contrast) / 50.0f;

                // reaction-diffusion
                float *tmp_cos = new float[num_strokes];
                float *tmp_sin = new float[num_strokes];
                for (int cycle = 0; cycle < reaction_diffusion_cycles; cycle++)
                {
                    for (int i = 0; i < num_strokes; i++)
                    {
                        const int n = stroke_offset(i);
                        strokes[i].rd_theta(neighbor_weight,
                                            refSAL[n], refG_Y[n],-refG_X[n],
                                            tmp_cos[i], tmp_sin[i]);
                    }
                    for (int i = 0; i < num_strokes; i++)
                    {
                        strokes[i].cos_theta = tmp_cos[i];
                        strokes[i].sin_theta = tmp_sin[i];
                    }
                }
                delete[] tmp_sin;
                delete[] tmp_cos;

                // update graph topology according to strokes' actual orientations
                compute_graph_topology();
            }

            // must redo all future steps
            ctvLEN = -1;
            ctvTHK = -1;
            ctvLIT = -1;
            ctvCHR = -1;
            ctvHUE = -1;
        }
    }

    void compute_graph_topology()
    {
        for (int i = 0; i < num_strokes; i++)
        {
            strokes[i].do_nn_clear();
            for (int j = 0; j < num_strokes; j++)
            {
                if (j != i)
                    strokes[i].do_nn_check(&strokes[j]);
            }
        }
    }

    void compute_size_u(const float min_u, const float max_u, const int contrast)
    {
        if (abs(minLEN - min_u) > 1e-6f || abs(maxLEN - max_u) > 1e-6f || ctvLEN != contrast)
        {
            minLEN = min_u;
            maxLEN = max_u;
            ctvLEN = contrast;

            // assign length values according to the reference map
            for (int i = 0; i < num_strokes; i++)
                strokes[i].u = minLEN + (maxLEN - minLEN) * refSIZ[stroke_offset(i)];

            if (contrast != 50)
            {
                const float neighbor_weight = (float)(50 - contrast) / 100.0f;

                // reaction-diffusion
                float *tmp_u = new float[num_strokes];
                for (int cycle = 0; cycle < reaction_diffusion_cycles; cycle++)
                {
                    for (int i = 0; i < num_strokes; i++)
                    {
                        const int n = stroke_offset(i);
                        tmp_u[i] = strokes[i].rd_u(neighbor_weight,
                                                   refSAL[n], minLEN + (maxLEN - minLEN) * refSIZ[n],
                                                   minLEN, maxLEN);
                    }
                    for (int i = 0; i < num_strokes; i++)
                        strokes[i].u = tmp_u[i];
                }
                delete[] tmp_u;
            }
        }
    }

    void compute_size_v(const float min_v, const float max_v, const int contrast)
    {
        if (abs(minTHK - min_v) > 1e-6f || abs(maxTHK - max_v) > 1e-6f || ctvTHK != contrast)
        {
            minTHK = min_v;
            maxTHK = max_v;
            ctvTHK = contrast;

            // assign thickness values according to the reference map
            for (int i = 0; i < num_strokes; i++)
                strokes[i].v = minTHK + (maxTHK - minTHK) * refSIZ[stroke_offset(i)];

            if (contrast != 50)
            {
                const float neighbor_weight = (float)(50 - contrast) / 100.0f;

                // reaction-diffusion
                float *tmp_v = new float[num_strokes];
                for (int cycle = 0; cycle < reaction_diffusion_cycles; cycle++)
                {
                    for (int i = 0; i < num_strokes; i++)
                    {
                        const int n = stroke_offset(i);
                        tmp_v[i] = strokes[i].rd_v(neighbor_weight,
                                                   refSAL[n], minTHK + (maxTHK - minTHK) * refSIZ[n],
                                                   minTHK, maxTHK);
                    }
                    for (int i = 0; i < num_strokes; i++)
                        strokes[i].v = tmp_v[i];
                }
                delete[] tmp_v;
            }
        }
    }

    void compute_color_L(const int contrast)
    {
        if (ctvLIT != contrast)
        {
            ctvLIT = contrast;

            // assign lightness values according to the source image
            for (int i = 0; i < num_strokes; i++)
                strokes[i].L = refLIT[stroke_offset(i)];

            if (contrast != 50)
            {
                const float neighbor_weight = (float)(50 - contrast) / 100.0f;
                const float prior_weight = 0.25f;

                // reaction-diffusion
                float *tmp_L = new float[num_strokes];
                for (int cycle = 0; cycle < reaction_diffusion_cycles; cycle++)
                {
                    for (int i = 0; i < num_strokes; i++)
                    {
                        tmp_L[i] = strokes[i].rd_L(neighbor_weight,
                                                   prior_weight, refLIT[stroke_offset(i)]);
                    }
                    for (int i = 0; i < num_strokes; i++)
                        strokes[i].L = tmp_L[i];
                }
                delete[] tmp_L;
            }
        }
    }

    void compute_color_C(const int contrast)
    {
        if (ctvCHR != contrast)
        {
            ctvCHR = contrast;

            // assign chroma values according to the source image
            for (int i = 0; i < num_strokes; i++)
                strokes[i].C = refCHR[stroke_offset(i)];

            if (contrast != 50)
            {
                const float neighbor_weight = (float)(50 - contrast) / 100.0f;
                const float prior_weight = 0.25f;

                // reaction-diffusion
                float *tmp_C = new float[num_strokes];
                for (int cycle = 0; cycle < reaction_diffusion_cycles; cycle++)
                {
                    for (int i = 0; i < num_strokes; i++)
                    {
                        tmp_C[i] = strokes[i].rd_C(neighbor_weight,
                                                   prior_weight, refCHR[stroke_offset(i)]);
                    }
                    for (int i = 0; i < num_strokes; i++)
                        strokes[i].C = tmp_C[i];
                }
                delete[] tmp_C;
            }
        }
    }

    void compute_color_H(const int contrast)
    {
        if (ctvHUE != contrast)
        {
            ctvHUE = contrast;

            // assign hue values according to the source image
            for (int i = 0; i < num_strokes; i++)
            {
                const int n = stroke_offset(i);
                strokes[i].cos_H = refC_H[n];
                strokes[i].sin_H = refS_H[n];
            }

            if (contrast != 50)
            {
                const float neighbor_weight = (float)(50 - contrast) / 50.0f;
                const float prior_weight = 0.25f;

                // reaction-diffusion
                float *tmp_cos = new float[num_strokes];
                float *tmp_sin = new float[num_strokes];
                for (int cycle = 0; cycle < reaction_diffusion_cycles; cycle++)
                {
                    for (int i = 0; i < num_strokes; i++)
                    {
                        const int n = stroke_offset(i);
                        strokes[i].rd_H(neighbor_weight,
                                        prior_weight, refC_H[n], refS_H[n],
                                        tmp_cos[i], tmp_sin[i]);
                    }
                    for (int i = 0; i < num_strokes; i++)
                    {
                        strokes[i].cos_H = tmp_cos[i];
                        strokes[i].sin_H = tmp_sin[i];
                    }
                }
                delete[] tmp_cos;
                delete[] tmp_sin;
            }
        }
    }

    inline int stroke_offset(int stroke_index)
    {
        return strokes[stroke_index].x + strokes[stroke_index].y * width;
    }
};

#endif // STROKEPATTERN_H
