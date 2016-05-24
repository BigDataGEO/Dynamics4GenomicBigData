{\rtf1\ansi\ansicpg1252\cocoartf1348\cocoasubrtf170
{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red13\green0\blue129;\red235\green236\blue237;\red36\green38\blue41;
\red37\green127\blue159;\red114\green121\blue129;\red104\green26\blue29;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720

\f0\fs26 \cf2 \cb3 \expnd0\expndtw0\kerning0
package\cf4 \expnd0\expndtw0\kerning0
 com.stackoverflow;\
\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 java.io.\cf5 \expnd0\expndtw0\kerning0
BufferedReader\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 java.io.\cf5 \expnd0\expndtw0\kerning0
InputStream\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 java.io.\cf5 \expnd0\expndtw0\kerning0
InputStreamReader\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 java.net.\cf5 \expnd0\expndtw0\kerning0
HttpURLConnection\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 java.net.URL;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 java.security.cert.X509Certificate;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 javax.net.ssl.\cf5 \expnd0\expndtw0\kerning0
HttpsURLConnection\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 javax.net.ssl.\cf5 \expnd0\expndtw0\kerning0
SSLContext\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 javax.net.ssl.\cf5 \expnd0\expndtw0\kerning0
SSLSession\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 javax.net.ssl.\cf5 \expnd0\expndtw0\kerning0
TrustManager\cf4 \expnd0\expndtw0\kerning0
;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 javax.net.ssl.X509TrustManager;\
\cf2 \expnd0\expndtw0\kerning0
import\cf4 \expnd0\expndtw0\kerning0
 javax.net.ssl.\cf5 \expnd0\expndtw0\kerning0
HostnameVerifier\cf4 \expnd0\expndtw0\kerning0
;\
\
\cf2 \expnd0\expndtw0\kerning0
public\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
class\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
Downloader\cf4 \expnd0\expndtw0\kerning0
 \{\
    \cf2 \expnd0\expndtw0\kerning0
public\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
static\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
 getData(\cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
 address) \cf2 \expnd0\expndtw0\kerning0
throws\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
Exception\cf4 \expnd0\expndtw0\kerning0
 \{\
        \cf6 \expnd0\expndtw0\kerning0
// Create a trust manager that does not validate certificate chains\cf4 \expnd0\expndtw0\kerning0
\
        \cf5 \expnd0\expndtw0\kerning0
TrustManager\cf4 \expnd0\expndtw0\kerning0
[] trustAllCerts = \cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
TrustManager\cf4 \expnd0\expndtw0\kerning0
[] \{\
            \cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 X509TrustManager() \{\
                \cf2 \expnd0\expndtw0\kerning0
public\cf4 \expnd0\expndtw0\kerning0
 java.security.cert.X509Certificate[] getAcceptedIssuers() \{\
                    \cf2 \expnd0\expndtw0\kerning0
return\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
null\cf4 \expnd0\expndtw0\kerning0
;\
                \}\
                \cf2 \expnd0\expndtw0\kerning0
public\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
void\cf4 \expnd0\expndtw0\kerning0
 checkClientTrusted(X509Certificate[] certs, \cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
 authType) \{\
                \}\
                \cf2 \expnd0\expndtw0\kerning0
public\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
void\cf4 \expnd0\expndtw0\kerning0
 checkServerTrusted(X509Certificate[] certs, \cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
 authType) \{\
                \}\
            \}\
        \};\
\
        \cf6 \expnd0\expndtw0\kerning0
// Create a host name verifier that always passes\cf4 \expnd0\expndtw0\kerning0
\
        \cf5 \expnd0\expndtw0\kerning0
HostnameVerifier\cf4 \expnd0\expndtw0\kerning0
 allHostsValid = \cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
HostnameVerifier\cf4 \expnd0\expndtw0\kerning0
() \{\
            \cf2 \expnd0\expndtw0\kerning0
public\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
boolean\cf4 \expnd0\expndtw0\kerning0
 verify(\cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
 hostname, \cf5 \expnd0\expndtw0\kerning0
SSLSession\cf4 \expnd0\expndtw0\kerning0
 session) \{\
                \cf2 \expnd0\expndtw0\kerning0
return\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
true\cf4 \expnd0\expndtw0\kerning0
;\
            \}\
        \};\
\
        \cf6 \expnd0\expndtw0\kerning0
// Install the all-trusting trust manager\cf4 \expnd0\expndtw0\kerning0
\
        \cf5 \expnd0\expndtw0\kerning0
SSLContext\cf4 \expnd0\expndtw0\kerning0
 sc = \cf5 \expnd0\expndtw0\kerning0
SSLContext\cf4 \expnd0\expndtw0\kerning0
.getInstance(\cf7 \expnd0\expndtw0\kerning0
"SSL"\cf4 \expnd0\expndtw0\kerning0
);\
        sc.init(\cf2 \expnd0\expndtw0\kerning0
null\cf4 \expnd0\expndtw0\kerning0
, trustAllCerts, \cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 java.security.\cf5 \expnd0\expndtw0\kerning0
SecureRandom\cf4 \expnd0\expndtw0\kerning0
());\
        \cf5 \expnd0\expndtw0\kerning0
HttpsURLConnection\cf4 \expnd0\expndtw0\kerning0
.setDefaultSSLSocketFactory(sc.getSocketFactory());\
\
        \cf6 \expnd0\expndtw0\kerning0
// Install the all-trusting host verifier\cf4 \expnd0\expndtw0\kerning0
\
        \cf5 \expnd0\expndtw0\kerning0
HttpsURLConnection\cf4 \expnd0\expndtw0\kerning0
.setDefaultHostnameVerifier(allHostsValid);\
\
        \cf6 \expnd0\expndtw0\kerning0
// open connection\cf4 \expnd0\expndtw0\kerning0
\
        URL page = \cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 URL(address);\
        \cf5 \expnd0\expndtw0\kerning0
HttpURLConnection\cf4 \expnd0\expndtw0\kerning0
 conn = (\cf5 \expnd0\expndtw0\kerning0
HttpURLConnection\cf4 \expnd0\expndtw0\kerning0
) page.openConnection();\
        \cf5 \expnd0\expndtw0\kerning0
BufferedReader\cf4 \expnd0\expndtw0\kerning0
 buff = \cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
BufferedReader\cf4 \expnd0\expndtw0\kerning0
(\cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
InputStreamReader\cf4 \expnd0\expndtw0\kerning0
(conn.getInputStream()));\
\
        \cf6 \expnd0\expndtw0\kerning0
// read text\cf4 \expnd0\expndtw0\kerning0
\
        \cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
 line;\
        \cf5 \expnd0\expndtw0\kerning0
StringBuffer\cf4 \expnd0\expndtw0\kerning0
 text = \cf2 \expnd0\expndtw0\kerning0
new\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
StringBuffer\cf4 \expnd0\expndtw0\kerning0
();\
        \cf2 \expnd0\expndtw0\kerning0
while\cf4 \expnd0\expndtw0\kerning0
 ( (line = buff.readLine()) != \cf2 \expnd0\expndtw0\kerning0
null\cf4 \expnd0\expndtw0\kerning0
 ) \{\
            \cf6 \expnd0\expndtw0\kerning0
//System.out.println(line);\cf4 \expnd0\expndtw0\kerning0
\
            text.append(line + \cf7 \expnd0\expndtw0\kerning0
"\\n"\cf4 \expnd0\expndtw0\kerning0
);\
        \}\
        buff.close();\
\
        \cf2 \expnd0\expndtw0\kerning0
return\cf4 \expnd0\expndtw0\kerning0
 text.toString();\
    \}\
\
    \cf2 \expnd0\expndtw0\kerning0
public\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
static\cf4 \expnd0\expndtw0\kerning0
 \cf2 \expnd0\expndtw0\kerning0
void\cf4 \expnd0\expndtw0\kerning0
 main(\cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
[] argv) \cf2 \expnd0\expndtw0\kerning0
throws\cf4 \expnd0\expndtw0\kerning0
 \cf5 \expnd0\expndtw0\kerning0
Exception\cf4 \expnd0\expndtw0\kerning0
 \{\
        \cf5 \expnd0\expndtw0\kerning0
String\cf4 \expnd0\expndtw0\kerning0
 str = getData(\cf7 \expnd0\expndtw0\kerning0
"https://expired.badssl.com/"\cf4 \expnd0\expndtw0\kerning0
);\
        \cf5 \expnd0\expndtw0\kerning0
System\cf4 \expnd0\expndtw0\kerning0
.out.println(str);\
    \}\
\}\
}