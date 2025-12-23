#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mesh.h"

//#define MAX_TOKENS 200
//#define MAX_STACK 200
/*
typedef enum {
    TOK_NUM,      // 숫자 (항상 양수 double)
    TOK_VAR,
    TOK_OP_BINARY,
    TOK_OP_UNARY, // 단항 + 또는 -
    TOK_FUNC,
    TOK_LPAREN,
    TOK_RPAREN
} TokenType;

typedef struct {
    TokenType type;
    union {
        double num;     // 항상 양수
        char var;
        char op;        // binary: + - * / ^, unary: + -
        char func[4];
    } value;
} Token;
*/

int precedence(char op) {
    if (op == '+' || op == '-') return 1;
    if (op == '*' || op == '/') return 2;
    if (op == '^') return 3;
    return 0;
}

int is_left_associative(char op) {
    return op != '^';
}

// 개선된 토큰화: 단항 +/- 를 반드시 별도 토큰으로 분리
int tokenize(const char* expr, Token* tokens) {
    int i = 0, idx = 0;
    while (expr[i]) {
        if (isspace(expr[i])) { i++; continue; }

        // 단항 + 또는 - 체크 (수식 시작, 연산자 후, (, 함수 후)
        char prev = (idx == 0) ? ' ' : (tokens[idx-1].type == TOK_OP_BINARY || 
                                        tokens[idx-1].type == TOK_OP_UNARY || 
                                        tokens[idx-1].type == TOK_LPAREN || 
                                        tokens[idx-1].type == TOK_FUNC) ? ' ' : '0';

        if ((expr[i] == '+' || expr[i] == '-') && prev == ' ') {
            tokens[idx].type = TOK_OP_UNARY;
            tokens[idx].value.op = expr[i];
            idx++;
            i++;
            // 바로 다음에 숫자가 없으면 에러 (예: "1 + -")
            if (!isdigit(expr[i]) && expr[i] != '.') {
                printf("Expected number after unary %c\n", expr[i-1]);
                exit(1);
            }
        }

        // 숫자 파싱 (항상 양수로, 지수 표기 지원)
        if (isdigit(expr[i]) || expr[i] == '.') {
            char* endptr;
            double val = strtod(&expr[i], &endptr);
            if (endptr == &expr[i]) {
                printf("Invalid number\n");
                exit(1);
            }
            tokens[idx].type = TOK_NUM;
            tokens[idx].value.num = val;  // 양수 값
            i = endptr - expr;
            idx++;
            continue;
        }

        if (expr[i] == 'x' || expr[i] == 'y') {
            tokens[idx].type = TOK_VAR;
            tokens[idx].value.var = expr[i++];
            idx++;
            continue;
        }

        if (expr[i] == '(') { tokens[idx].type = TOK_LPAREN; i++; idx++; continue; }
        if (expr[i] == ')') { tokens[idx].type = TOK_RPAREN; i++; idx++; continue; }

        if (strchr("+-*/^", expr[i])) {
            tokens[idx].type = TOK_OP_BINARY;
            tokens[idx].value.op = expr[i++];
            idx++;
            continue;
        }

        if (isalpha(expr[i])) {
            char func[4] = {0};
            int j = 0;
            while (isalpha(expr[i]) && j < 3) func[j++] = expr[i++];
            if (strcmp(func, "sin") == 0 || 
                strcmp(func, "cos") == 0 || 
                strcmp(func, "tan") == 0 || 
                strcmp(func, "exp") == 0) {
                tokens[idx].type = TOK_FUNC;
                strcpy(tokens[idx].value.func, func);
                idx++;
                continue;
            }
            printf("Unknown function: %s\n", func);
            exit(1);
        }

        printf("Invalid character: %c\n", expr[i]);
        exit(1);
    }
    return idx;
}

// Shunting-yard 및 evaluate_rpn은 이전과 거의 동일 (단항 우선순위 높게 처리)
void infix_to_rpn(Token* infix, int infix_len, Token* rpn, int* rpn_len) {
    Token op_stack[MAX_STACK];
    int op_top = -1;
    *rpn_len = 0;

    for (int i = 0; i < infix_len; i++) {
        Token t = infix[i];
        if (t.type == TOK_NUM || t.type == TOK_VAR) {
            rpn[(*rpn_len)++] = t;
        } else if (t.type == TOK_FUNC || t.type == TOK_OP_UNARY) {
            op_stack[++op_top] = t;  // 함수와 단항은 최고 우선순위
        } else if (t.type == TOK_LPAREN) {
            op_stack[++op_top] = t;
        } else if (t.type == TOK_RPAREN) {
            while (op_top >= 0 && op_stack[op_top].type != TOK_LPAREN) {
                rpn[(*rpn_len)++] = op_stack[op_top--];
            }
            if (op_top >= 0) op_top--;
            if (op_top >= 0 && op_stack[op_top].type == TOK_FUNC) {
                rpn[(*rpn_len)++] = op_stack[op_top--];
            }
        } else if (t.type == TOK_OP_BINARY) {
            while (op_top >= 0 && 
                   (op_stack[op_top].type == TOK_OP_BINARY || op_stack[op_top].type == TOK_OP_UNARY) &&
                   precedence(op_stack[op_top].value.op) >= precedence(t.value.op) &&
                   (precedence(op_stack[op_top].value.op) > precedence(t.value.op) ||
                    is_left_associative(t.value.op))) {
                rpn[(*rpn_len)++] = op_stack[op_top--];
            }
            op_stack[++op_top] = t;
        }
    }
    while (op_top >= 0) {
        if (op_stack[op_top].type == TOK_LPAREN) { printf("Mismatched parentheses\n"); exit(1); }
        rpn[(*rpn_len)++] = op_stack[op_top--];
    }
}

double evaluate_rpn(Token* rpn, int rpn_len, double x_val, double y_val) {
    double stack[MAX_STACK];
    int top = -1;

    for (int i = 0; i < rpn_len; i++) {
        Token t = rpn[i];
        if (t.type == TOK_NUM) {
            stack[++top] = t.value.num;
        } else if (t.type == TOK_VAR) {
            stack[++top] = (t.value.var == 'x') ? x_val : y_val;
        } else if (t.type == TOK_FUNC) {
            double a = stack[top--];
            if (strcmp(t.value.func, "sin") == 0) stack[++top] = sin(a);
            else if (strcmp(t.value.func, "cos") == 0) stack[++top] = cos(a);
            else if (strcmp(t.value.func, "tan") == 0) stack[++top] = tan(a);
            else stack[++top] = exp(a);
        } else if (t.type == TOK_OP_UNARY) {
            double a = stack[top--];
            stack[++top] = (t.value.op == '-') ? -a : a;
        } else if (t.type == TOK_OP_BINARY) {
            double b = stack[top--];
            double a = stack[top--];
            switch (t.value.op) {
                case '+': stack[++top] = a + b; break;
                case '-': stack[++top] = a - b; break;
                case '*': stack[++top] = a * b; break;
                case '/': if (b == 0.0) { printf("Division by zero\n"); exit(1); } stack[++top] = a / b; break;
                case '^': stack[++top] = pow(a, b); break;
            }
        }
    }
    return stack[top];
}

/*
int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: %s \"expression\" x_value y_value\n", argv[0]);
        return 1;
    }

    const char* expr = argv[1];
    double x = strtod(argv[2], NULL);
    double y = strtod(argv[3], NULL);

    Token infix[MAX_TOKENS];
    int infix_len = tokenize(expr, infix);

    Token rpn[MAX_TOKENS];
    int rpn_len;
    infix_to_rpn(infix, infix_len, rpn, &rpn_len);

    double result = evaluate_rpn(rpn, rpn_len, x, y);
    printf("Result: %.10g\n", result);

    return 0;
}
*/
