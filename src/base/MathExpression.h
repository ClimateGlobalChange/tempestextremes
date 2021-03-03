///////////////////////////////////////////////////////////////////////////////
///
///	\file    MathExpression.h
///	\author  t-mat, modified by Paul Ullrich
///	\version March 2nd, 2021
///
///	<remarks>
///		Based on http://en.wikipedia.org/wiki/Shunting-yard_algorithm
///		Code by https://ideone.com/DYX5CW
///	</remarks>

#ifndef _MATHEXPRESSION_H_
#define _MATHEXPRESSION_H_

#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <deque>
#include <cstdio>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

class MathExpression {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	MathExpression(
		const std::string & str
	) :
		m_fPostfix(false)
	{
		Tokenize(str);
	}

public:
	///	<summary>
	///		A class representing a token in an expression.
	///	</summary>
	class Token {
	public:
		enum class Type {
			Unknown,
			Number,
			Variable,
			Operator,
			LeftParen,
			RightParen,
		};

		Token(
			Type t,
			const std::string & s,
			int prec = -1,
			bool ra = false
		) :
			type(t),
			str(s),
			precedence(prec),
			rightAssociative(ra)
		{ }

		Type type;
		std::string str;
		int precedence;
		bool rightAssociative;
	};

protected:
	///	<summary>
	///		Debug Report.
	///	</summary>
	template<class T, class U>
	void debugReport(
		const Token& token,
		const T& queue,
		const U& stack,
		const std::string& comment = ""
	) {
		std::ostringstream ossQueue;
		for(const auto& t : queue) {
			ossQueue << " " << t.str;
		}
	
		std::ostringstream ossStack;
		for(const auto& t : stack) {
			ossStack << " " << t.str;
		}
	
		printf("|%-3s|%-32s|%10s| %s\n"
			   , token.str.c_str()
			   , ossQueue.str().c_str()
			   , ossStack.str().c_str()
			   , comment.c_str()
		);
	}

	///	<summary>
	///		Tokenize an input expression.
	///	</summary>
	void Tokenize(const std::string& expr) {
	
		for(const auto* p = expr.c_str(); *p; ++p) {
			if ((*p) == ' ') {
				continue;

			} else if (isdigit(*p)) {
				const auto* b = p;
				for(; isdigit(*p); ++p) {
					;
				}
				const auto s = std::string(b, p);
				m_tokens.push_back(Token(Token::Type::Number, s));
				--p;
	
			} else if (isalpha(*p)) {
				const auto* b = p;
				for(; isalpha(*p) || isdigit(*p); ++p) {
					;
				}
				const auto s = std::string(b, p);
				m_tokens.push_back(Token(Token::Type::Variable, s));
				--p;
	
			} else {
				Token::Type t = Token::Type::Unknown;
				int pr = -1;
				bool ra = false;
				switch(*p) {
				default:                                 break;
				case '(':   t = Token::Type::LeftParen;	 break;
				case ')':   t = Token::Type::RightParen; break;
				case '^':   t = Token::Type::Operator;    pr = 4; ra = true;  break;
				case '*':   t = Token::Type::Operator;    pr = 3; break;
				case '/':   t = Token::Type::Operator;    pr = 3; break;
				case '+':   t = Token::Type::Operator;    pr = 2; break;
				case '-':   t = Token::Type::Operator;    pr = 2; break;
				}
				m_tokens.push_back(Token(t, std::string(1, *p), pr, ra));
			}
		}
	}

public:
	///	<summary>
	///		Perform the shunting yard algorithm to convert an infix string to postfix.
	///	</summary>
	void ShuntingYard() {
		if (m_fPostfix) {
			_EXCEPTIONT("MathExpression already in postfix form");
		}

		std::deque<Token> queue;
		std::vector<Token> stack;
	
		// While there are tokens to be read:
		for(auto token : m_tokens) {
			// Read a token
			switch(token.type) {
			case Token::Type::Number:
			case Token::Type::Variable:
				// If the token is a number, then add it to the output queue
				queue.push_back(token);
				break;

			case Token::Type::Operator:
				{
					// If the token is operator, o1, then:
					const auto o1 = token;
	
					// while there is an operator token,
					while(!stack.empty()) {
						// o2, at the top of stack, and
						const auto o2 = stack.back();
	
						// either o1 is left-associative and its precedence is
						// *less than or equal* to that of o2,
						// or o1 if right associative, and has precedence
						// *less than* that of o2,
						if((! o1.rightAssociative && o1.precedence <= o2.precedence) ||
						   (  o1.rightAssociative && o1.precedence <  o2.precedence)
						) {
							// then pop o2 off the stack,
							stack.pop_back();
							// onto the output queue;
							queue.push_back(o2);
	
							continue;
						}
	
						// @@ otherwise, exit.
						break;
					}
	
					// push o1 onto the stack.
					stack.push_back(o1);
				}
				break;
	
			case Token::Type::LeftParen:
				// If token is left parenthesis, then push it onto the stack
				stack.push_back(token);
				break;
	
			case Token::Type::RightParen:
				// If token is right parenthesis:
				{
					bool match = false;
	
					// Until the token at the top of the stack
					// is a left parenthesis,
					while(! stack.empty() && stack.back().type != Token::Type::LeftParen) {
						// pop operators off the stack
						// onto the output queue.
						queue.push_back(stack.back());
						stack.pop_back();
						match = true;
					}
	
					// Pop the left parenthesis from the stack,
					// but not onto the output queue.
					stack.pop_back();
	
					if(!match && stack.empty()) {
						// If the stack runs out without finding a left parenthesis,
						// then there are mismatched parentheses.
						_EXCEPTION1("RightParen error (%s)", token.str.c_str());
					}
				}
				break;
	
			default:
				_EXCEPTION1("error (%s)", token.str.c_str());
			}
	
			//debugReport(token, queue, stack);
		}
	
		// When there are no more tokens to read:
		//   While there are still operator tokens in the stack:
		while(! stack.empty()) {
			// If the operator token on the top of the stack is a parenthesis,
			// then there are mismatched parentheses.
			if(stack.back().type == Token::Type::LeftParen) {
				_EXCEPTIONT("Mismatched parentheses error");
			}
	
			// Pop the operator onto the output queue.
			queue.push_back(std::move(stack.back()));
			stack.pop_back();
		}
	
		//debugReport(Token(Token::Type::Unknown, "End"), queue, stack);
	
		m_fPostfix = true;
		m_tokens = queue;
	}
	
///////////////////////////////////////////////////////////////////////////////

	//void Evaluate(
	//	std::map<std::string, double>
	/*
		for(const auto& expr : expressions) {
			printf("expr = %s\n", expr.c_str());
			printf("|%-3s|%-32s|%-10s|\n", "Tkn", "Queue", "Stack");
	
			const auto tokens = exprToTokens(expr);
			auto queue = shuntingYard(tokens);
			std::vector<int> stack;
	
			printf("\nCalculation\n");
			printf("|%-3s|%-32s|%-10s|\n", "Tkn", "Queue", "Stack");
	
			while(! queue.empty()) {
				std::string op;
	
				const auto token = queue.front();
				queue.pop_front();
				switch(token.type) {
				case Token::Type::Number:
					stack.push_back(std::stoi(token.str));
					op = "Push " + token.str;
					break;
	
				case Token::Type::Operator:
					{
						const auto rhs = stack.back();
						stack.pop_back();
						const auto lhs = stack.back();
						stack.pop_back();
	
						switch(token.str[0]) {
						default:
							printf("Operator error [%s]\n", token.str.c_str());
							exit(0);
							break;
						case '^':
							stack.push_back(static_cast<int>(pow(lhs, rhs)));
							break;
						case '*':
							stack.push_back(lhs * rhs);
							break;
						case '/':
							stack.push_back(lhs / rhs);
							break;
						case '+':
							stack.push_back(lhs + rhs);
							break;
						case '-':
							stack.push_back(lhs - rhs);
							break;
						}
						op = "Push " + std::to_string(lhs) + " " + token.str + " " + std::to_string(rhs);
					}
					break;
	
				default:
					printf("Token error\n");
					exit(0);
				}
				debugReport(token, queue, stack, op);
			}
			printf("\n  result = %d\n\n", stack.back());
		}
	}
*/

public:
	///	<summary>
	///		Number of tokens in the expression.
	///	</summary>
	size_t size() const {
		return m_tokens.size();
	}

	///	<summary>
	///		Access element at specific location.
	///	</summary>
	Token & operator[](size_t s) {
		return m_tokens[s];
	}

	///	<summary>
	///		Access element at specific location.
	///	</summary>
	const Token & operator[](size_t s) const {
		return m_tokens[s];
	}

public:
	///	<summary>
	///		Boolean indicating the tokens are in postfix notation.
	///	</summary>
	bool m_fPostfix;

	///	<summary>
	///		Tokenized expression.
	///	</summary>
	std::deque<Token> m_tokens;

///////////////////////////////////////////////////////////////////////////////

};

#endif // _MATHEXPRESSION_H_

