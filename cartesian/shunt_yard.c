#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define MAXOPSTACK 64
#define MAXNUMSTACK 64
#define MAXRPN 64  // Max tokens in RPN

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

typedef enum {
    TOK_NUM,
    TOK_VAR,
    TOK_OP
} TokenType;

// Token structure
typedef struct {
    TokenType type;
    double num;  // For TOK_NUM
    struct op_s *op;  // For TOK_OP
} Token;


// 연산자/함수 정의
enum {ASSOC_NONE=0, ASSOC_LEFT, ASSOC_RIGHT};

// Binary operation functions
static double op_add(double a, double b) { return a + b; }
static double op_sub(double a, double b) { return a - b; }
static double op_mul(double a, double b) { return a * b; }
static double op_div(double a, double b) { return a / b; }

// Unary wrappers (ignore second arg to match signature)
static double wrap_sin(double a, double unused) { (void)unused; return sin(a); }
static double wrap_cos(double a, double unused) { (void)unused; return cos(a); }
static double wrap_exp(double a, double unused) { (void)unused; return exp(a); }

struct op_s {
    const char *name;  // 연산자 또는 함수 이름
    int prec;
    int assoc;
    int unary;
    double (*eval)(double a1, double a2);
} ops[] = {
    {"_", 10, ASSOC_RIGHT, 1, NULL},  // unary minus (eval 함수 별도 처리)
    {"^", 9, ASSOC_RIGHT, 0, pow},
    {"*", 8, ASSOC_LEFT, 0, op_mul},
    {"/", 8, ASSOC_LEFT, 0, op_div},
    {"+", 5, ASSOC_LEFT, 0, op_add},
    {"-", 5, ASSOC_LEFT, 0, op_sub},
    {"sin", 11, ASSOC_RIGHT, 1, wrap_sin},  // sin 함수 추가 (high prec, unary)
    {"cos", 11, ASSOC_RIGHT, 1, wrap_cos},  // cos 함수 추가 (high prec, unary)
    {"exp", 11, ASSOC_RIGHT, 1, wrap_exp},  // exp 함수 추가 (high prec, unary)
    {"(", 0, ASSOC_NONE, 0, NULL},
    {")", 0, ASSOC_NONE, 0, NULL},
    {NULL}  // 종료
};

struct op_s *getop(const char *s) {
    for (int i = 0; ops[i].name; ++i) {
        if (strcmp(ops[i].name, s) == 0) return &ops[i];
    }
    return NULL;
}

// Shunting-yard to RPN
int shunting_yard(const char *input, Token *rpn) {
    struct op_s *opstack[MAXOPSTACK];
    int nopstack = 0;
    int rpn_idx = 0;
    const char *p = input;
    char buf[32];
    int i;

    while (*p) {
        while (isspace(*p)) ++p;
        if (!*p) break;

        if (isdigit(*p) || *p == '.') {  // 숫자 (지수 표기 포함)
            i = 0;
				//정수/소수 부분
				if (isdigit(*p)) {
            	while (isdigit(*p)) buf[i++] = *p++;
            }
            if (*p == '.') {
            	buf[i++] = *p++;
            	while (isdigit(*p)) buf[i++] = *p++;
            }
				// 지수 부분 (e/E, 부호, 숫자)
            if (*p == 'e' || *p == 'E') {
            	buf[i++] = *p++;
            	if (*p == '+' || *p == '-') buf[i++] = *p++;
            	while (isdigit(*p)) buf[i++] = *p++;
            }

            //do buf[i++] = *p++; while (isdigit(*p) || *p == '.');
            buf[i] = '\0';
            rpn[rpn_idx].type = TOK_NUM;
            rpn[rpn_idx++].num = atof(buf);
        } else if (isalpha(*p)) {  // 함수 또는 변수
            i = 0;
            do buf[i++] = *p++; while (isalpha(*p));
            buf[i] = '\0';
            if (strcmp(buf, "x") == 0) {  // 변수 x
            	rpn[rpn_idx].type = TOK_VAR;
            	rpn_idx++;
            } else if (strcmp(buf, "pi") == 0) {  // 상수 pi
            	rpn[rpn_idx].type = TOK_NUM;
            	rpn[rpn_idx++].num = M_PI;
				} else {  // 함수
                struct op_s *op = getop(buf);
                if (!op) { fprintf(stderr, "Unknown function: %s\n", buf); exit(1); }
                // Shunt op
                while (nopstack && op->prec <= opstack[nopstack - 1]->prec &&
                       (op->assoc != ASSOC_RIGHT || op->prec < opstack[nopstack - 1]->prec)) {
                    rpn[rpn_idx].type = TOK_OP;
                    rpn[rpn_idx++].op = opstack[--nopstack];
                }
                opstack[nopstack++] = op;
            }
        } else {  // 연산자
            buf[0] = *p++;
            buf[1] = 0;
            if (buf[0] == '-' && (p == input + 1 || (!isdigit(*(p-2)) && *(p-2) != ')' && *(p-2) != 'x'))) {  // unary minus
                buf[0] = '_';
            }
            struct op_s *op = getop(buf);
            if (!op) { fprintf(stderr, "Unknown operator: %s\n", buf); exit(1); }

            if (strcmp(op->name, "(") == 0) {
                opstack[nopstack++] = op;
            } else if (strcmp(op->name, ")") == 0) {
                while (nopstack && strcmp(opstack[nopstack - 1]->name, "(") != 0) {
                    rpn[rpn_idx].type = TOK_OP;
                    rpn[rpn_idx++].op = opstack[--nopstack];
                }
                if (nopstack == 0 || strcmp(opstack[--nopstack]->name, "(") != 0) {
                    fprintf(stderr, "Mismatched parentheses\n");
                    exit(1);
                }
            } else {
                while (nopstack && op->prec <= opstack[nopstack - 1]->prec &&
                       (op->assoc != ASSOC_RIGHT || op->prec < opstack[nopstack - 1]->prec)) {
                    rpn[rpn_idx].type = TOK_OP;
                    rpn[rpn_idx++].op = opstack[--nopstack];
                }
                opstack[nopstack++] = op;
            }
        }
    }

    // Flush op stack to RPN
    while (nopstack) {
        rpn[rpn_idx].type = TOK_OP;
        rpn[rpn_idx++].op = opstack[--nopstack];
    }

    return rpn_idx;  // Number of tokens in RPN
}

// Evaluate RPN with given x
double evaluate_rpn(Token *rpn, int rpn_size, double x) {
    double numstack[MAXNUMSTACK];
    int nnumstack = 0;

    for (int i = 0; i < rpn_size; i++) {
        Token t = rpn[i];
        if (t.type == TOK_NUM) {
            numstack[nnumstack++] = t.num;
        } else if (t.type == TOK_VAR) {
            numstack[nnumstack++] = x;
        } else if (t.type == TOK_OP) {
            struct op_s *op = t.op;
            double n1 = numstack[--nnumstack];
            if (op->unary) {
                if (strcmp(op->name, "_") == 0) {
                    numstack[nnumstack++] = -n1;
                } else {
                    numstack[nnumstack++] = op->eval(n1, 0);
                }
            } else {
                double n2 = numstack[--nnumstack];
                numstack[nnumstack++] = op->eval(n2, n1);
            }
        }
    }

    return numstack[--nnumstack];
}

/*
int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("사용법: %s <표현식1> [<표현식2>] <x값>\n", argv[0]);
        printf("예시: %s \"sin(3*x)+exp(x)\" 1.0\n", argv[0]);
        printf("또는: %s \"sin(3*x)+exp(x)\" \"2*x^2+3\" 1.0\n", argv[0]);
        return 1;
    }

    int num_expr = argc - 2;  // 표현식 개수 (마지막 argv는 x)
    double x = atof(argv[argc - 1]);

    printf("num_exp=%d, argc=%d\n",num_expr,argc);

    for (int expr_idx = 1; expr_idx <= num_expr; expr_idx++) {
        const char *input = argv[expr_idx];
        Token rpn[MAXRPN];
        int rpn_size = shunting_yard(input, rpn);
        
        printf("expr_idx=%d, input=%s, rpn_size=%d\n",expr_idx,input,rpn_size);
        double result = evaluate_rpn(rpn, rpn_size, x);
        printf("표현식 \"%s\" 의 결과 (x=%.2f): %f\n", input, x, result);
    }

    return 0;
}
*/
