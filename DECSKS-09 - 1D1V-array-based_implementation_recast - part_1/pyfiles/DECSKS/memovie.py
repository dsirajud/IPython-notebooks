from numpy import *


# if using the top shape
A1 = (120 + 10.)*40
A2 = (20)*(20 + 120 + 20.)

ytilde1 = ( (120 + 10)/2. - 10)
ytilde2 = -10 - 20./2

ybar = (ytilde1*A1 + ytilde2*A2) / (A1 + A2)

Q_top = ybar*(A1 + A2)

# if using bottom two shapes

A_right = 40*20.
ybar = -10 - 20 - 40/2.

Q_bot_right = ybar*A_right

Q_bot = 2*Q_bot_right # Q_bot_left = Q_bot_right



Q1_top = (130./2 - 10)*(130*40)
Q2_top = (-10 - 20./2)*(20)*(120 + 20 + 20)


######### find the neutrarl axis location

A_bot1 = 40*20.
A_bot2 = A_bot1

A_top1 = 40*(120 + 10.)
A_top2 = 20*(20 + 120 + 20.)

ytilde_bot1 = 40./2
ytilde_bot2 = ytilde_bot1

ytilde_top2 = 40 + 20./2
ytilde_top1 = 40 + 20 + 130./2

ybar = (ytilde_top1 * A_top1 + ytilde_top2*A_top2 +
        ytilde_bot1 * A_bot1 + ytilde_bot2 * A_bot2 ) / \
        (A_bot1 + A_bot2 + A_top1 + A_top2)

        # ybar = 84.2 mm from the bottom

# Q calc, top

ytilde1 = 40 + 20 + 130./2 - ybar
ytilde2 = 40./2 - ybar

A1 = A_top1
A2 = A_top2

A_prime = A1 + A2

ybar_prime = (ytilde1 * A1 + ytilde2 * A2) / (A1 + A2)

Q_top = ybar_prime * A_prime

# Q aclc, both

yprime_right = 40./2 - ybar
A_right = 40 * 20.

Q_bot = yprime_right * A_right * 2

