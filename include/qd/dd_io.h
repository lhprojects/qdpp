#ifndef _QD_DD_IO_H
#define _QD_DD_IO_H

#include "dd_real.h"
#include "dd_math.h"

std::ostream& operator<<(std::ostream& s, const dd_real& a);
std::istream& operator>>(std::istream& s, dd_real& a);

inline dd_real ldexp_(const dd_real& a, int exp)
{
    return dd_real(std::ldexp(a.x[0], exp), std::ldexp(a.x[1], exp));
}

/* Outputs the double-double number dd. */
inline std::ostream& operator<<(std::ostream& os, const dd_real& dd)
{
    bool showpos = (os.flags() & std::ios_base::showpos) != 0;
    bool uppercase = (os.flags() & std::ios_base::uppercase) != 0;
    return os << dd.to_string((int)/*HL: */os.precision(), (int)/*HL: */os.width(), os.flags(),
        showpos, uppercase, os.fill());
}

/* Reads in the double-double number a. */
inline std::istream& operator>>(std::istream& s, dd_real& a)
{
    char str[255];
    s >> str;
    dd_real::read(str, a);
    return s;
}

inline void dd_real::to_digits(char* s, int& expn, int precision) const
{
    int D = precision + 1;  /* number of digits to compute */

    dd_real r = abs(*this);
    int e;  /* exponent */
    int i, d;

    if (x[0] == 0.0) {
        /* this == 0.0 */
        expn = 0;
        for (i = 0; i < precision; i++) s[i] = '0';
        return;
    }

    /* First determine the (approximate) exponent. */
    e = to_int(std::floor(std::log10(std::abs(x[0]))));

    if (e < -300) {
        r *= dd_real(10.0) ^ 300;
        r /= dd_real(10.0) ^ (e + 300);
    } else if (e > 300) {
        r = ldexp_(r, -53);
        r /= dd_real(10.0) ^ e;
        r = ldexp_(r, 53);
    } else {
        r /= dd_real(10.0) ^ e;
    }

    /* Fix exponent if we are off by one */
    if (r >= 10.0) {
        r /= 10.0;
        e++;
    } else if (r < 1.0) {
        r *= 10.0;
        e--;
    }

    if (r >= 10.0 || r < 1.0) {
        qd::error("(dd_real::to_digits): can't compute exponent.");
        return;
    }

    /* Extract the digits */
    for (i = 0; i < D; i++) {
        d = static_cast<int>(r.x[0]);
        r -= d;
        r *= 10.0;

        s[i] = static_cast<char>(d + '0');
    }

    /* Fix out of range digits. */
    for (i = D - 1; i > 0; i--) {
        if (s[i] < '0') {
            s[i - 1]--;
            s[i] += 10;
        } else if (s[i] > '9') {
            s[i - 1]++;
            s[i] -= 10;
        }
    }

    if (s[0] <= '0') {
        qd::error("(dd_real::to_digits): non-positive leading digit.");
        return;
    }

    /* Round, handle carry */
    if (s[D - 1] >= '5') {
        s[D - 2]++;

        i = D - 2;
        while (i > 0 && s[i] > '9') {
            s[i] -= 10;
            s[--i]++;
        }
    }

    /* If first digit is 10, shift everything. */
    if (s[0] > '9') {
        e++;
        for (i = precision; i >= 2; i--) s[i] = s[i - 1];
        s[0] = '1';
        s[1] = '0';
    }

    s[precision] = 0;
    expn = e;
}

/* Writes the double-double number into the character array s of length len.
   The integer d specifies how many significant digits to write.
   The string s must be able to hold at least (d+8) characters.
   showpos indicates whether to use the + sign, and uppercase indicates
   whether the E or e is to be used for the exponent. */
inline void dd_real::write(char* s, int len, int precision,
    bool showpos, bool uppercase) const
{
    std::string str = to_string(precision, 0, std::ios_base::scientific, showpos, uppercase);

    // MSVC don't like user to use strncpy
    size_t ln = (std::min)(str.length(), (size_t)len - 1);
    memcpy(s, str.data(), ln);
    s[len - 1] = '\0';
}


inline void round_string(char* s, int precision, int* offset)
{
    /*
     Input string must be all digits or errors will occur.
     */

    int i;
    int D = precision;

    /* Round, handle carry */
    if (s[D - 1] >= '5') {
        s[D - 2]++;

        i = D - 2;
        while (i > 0 && s[i] > '9') {
            s[i] -= 10;
            s[--i]++;
        }
    }

    /* If first digit is 10, shift everything. */
    if (s[0] > '9') {
        // e++; // don't modify exponent here
        for (i = precision; i >= 2; i--) s[i] = s[i - 1];
        s[0] = '1';
        s[1] = '0';

        (*offset)++; // now offset needs to be increased by one
        precision++;
    }

    s[precision] = 0; // add terminator for array
}

inline std::string dd_real::to_string(int precision, int width, std::ios_base::fmtflags fmt,
    bool showpos, bool uppercase, char fill) const
{
    std::string s;
    bool fixed = (fmt & std::ios_base::fixed) != 0;
    bool sgn = true;
    int i, e = 0;

    if (isnan()) {
        s = uppercase ? "NAN" : "nan";
        sgn = false;
    } else {
        if (*this < 0.0)
            s += '-';
        else if (showpos)
            s += '+';
        else
            sgn = false;

        if (isinf()) {
            s += uppercase ? "INF" : "inf";
        } else if (*this == 0.0) {
            /* Zero case */
            s += '0';
            if (precision > 0) {
                s += '.';
                s.append(precision, '0');
            }
        } else {
            /* Non-zero case */
            int off = (fixed ? (1 + to_int(floor(log10(abs(*this))))) : 1);
            int d = precision + off;

            int d_with_extra = d;
            if (fixed)
                d_with_extra = (std::max)(60, d); // longer than the max accuracy for DD

            // highly special case - fixed mode, precision is zero, abs(*this) < 1.0
            // without this trap a number like 0.9 printed fixed with 0 precision prints as 0
            // should be rounded to 1.
            if (fixed && (precision == 0) && (abs(*this) < 1.0)) {
                if (abs(*this) >= 0.5)
                    s += '1';
                else
                    s += '0';

                return s;
            }

            // handle near zero to working precision (but not exactly zero)
            if (fixed && d <= 0) {
                s += '0';
                if (precision > 0) {
                    s += '.';
                    s.append(precision, '0');
                }
            } else { // default

                char* t; //  = new char[d+1];
                int j;

                if (fixed) {
                    t = new char[(size_t)d_with_extra + 1]; //HL: safe new
                    to_digits(t, e, d_with_extra);
                } else {
                    t = new char[(size_t)d + 1];   //HL: safe new
                    to_digits(t, e, d);
                }

                if (fixed) {
                    // fix the string if it's been computed incorrectly
                    // round here in the decimal string if required
                    round_string(t, d + 1, &off);

                    if (off > 0) {
                        for (i = 0; i < off; i++) s += t[i];
                        if (precision > 0) {
                            s += '.';
                            for (j = 0; j < precision; j++, i++) s += t[i];
                        }
                    } else {
                        s += "0.";
                        if (off < 0) s.append(-off, '0');
                        for (i = 0; i < d; i++) s += t[i];
                    }
                } else {
                    s += t[0];
                    if (precision > 0) s += '.';

                    for (i = 1; i <= precision; i++)
                        s += t[i];

                }
                delete[] t;
            }
        }

        // trap for improper offset with large values
        // without this trap, output of values of the for 10^j - 1 fail for j > 28
        // and are output with the point in the wrong place, leading to a dramatically off value
        if (fixed && (precision > 0)) {
            // make sure that the value isn't dramatically larger
            double from_string = atof(s.c_str());

            // if this ratio is large, then we've got problems
            if (fabs(from_string / this->x[0]) > 3.0) {

                //HL: unreferenced local variable   int point_position;
                //HL: unreferenced local variable   char temp;

                // loop on the string, find the point, move it up one
                // don't act on the first character
                for (i = 1; i < (int)/*HL: */s.length(); i++) {
                    if (s[i] == '.') {
                        s[i] = s[i - 1];
                        s[i - 1] = '.';
                        break;
                    }
                }

                from_string = atof(s.c_str());
                // if this ratio is large, then the string has not been fixed
                if (fabs(from_string / this->x[0]) > 3.0) {
                    qd::error("Re-rounding unsuccessful in large number fixed point trap.");
                }
            }
        }


        if (!fixed && !isinf()) {
            /* Fill in exponent part */
            s += uppercase ? 'E' : 'e';
            append_expn(s, e);
        }
    }

    /* Fill in the blanks */
    int len = (int)s.length(); //HL:suppress warnings
    if (len < width) {
        int delta = width - len;
        if (fmt & std::ios_base::internal) {
            if (sgn)
                s.insert(static_cast<std::string::size_type>(1), delta, fill);
            else
                s.insert(static_cast<std::string::size_type>(0), delta, fill);
        } else if (fmt & std::ios_base::left) {
            s.append(delta, fill);
        } else {
            s.insert(static_cast<std::string::size_type>(0), delta, fill);
        }
    }

    return s;
}


/* Debugging routines */
inline void dd_real::dump(const std::string& name, std::ostream& os) const
{
    std::ios_base::fmtflags old_flags = os.flags();
    std::streamsize old_prec = os.precision(19);
    os << std::scientific;

    if (name.length() > 0) os << name << " = ";
    os << "[ " << std::setw(27) << x[0] << ", " << std::setw(27) << x[1] << " ]" << std::endl;

    os.precision(old_prec);
    os.flags(old_flags);
}

inline void dd_real::dump_bits(const std::string& name, std::ostream& os) const
{
    std::string::size_type len = name.length();
    if (len > 0) {
        os << name << " = ";
        len += 3;
    }
    os << "[ ";
    len += 2;
    print_double_info(os, x[0]);
    os << std::endl;
    for (std::string::size_type i = 0; i < len; i++) os << ' ';
    print_double_info(os, x[1]);
    os << " ]" << std::endl;
}
#endif
