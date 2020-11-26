################################################################################
#         Métodos Numéricos para Equações Diferenciais II -- Trabalho 2        #
#                           Ariel Nogueira Kovaljski                           #
################################################################################

import matplotlib.pyplot as plt


x0, y0 = [], []
x1, y1 = [], []
x2, y2 = [], []
parameters = {}

def main():
    # Abre arquivos para leitura
    # Arquivo Upwind
    try:
        f0 = open('./results/results0.txt', 'r')
    except OSError as err:
        print("Erro ao abrir o arquivo 'results0.txt':", err)

    # Arquivo Superbee
    try:
        f1 = open('./results/results1.txt', 'r')
    except OSError as err:
        print("Erro ao abrir o arquivo 'results1.txt':", err)

    # Arquivo Van Albada
    try:
        f2 = open('./results/results2.txt', 'r')
    except OSError as err:
        print("Erro ao abrir o arquivo 'results2.txt':", err)


    # Extrai informações dos arquivos
    data_extract(f0, x0, y0)
    data_extract(f1, x1, y1)
    data_extract(f2, x2, y2)

    # Plota os gráficos de cada método
    plot_graph(x0, y0, 0)
    plot_graph(x1, y1, 1)
    plot_graph(x2, y2, 2)

    # Plota todos os gráficos na mesma janela
    plot_all_graphs()


def data_extract(f, x, y):
    for line_number, line in enumerate(f):
        if (line_number + 1) < 6:
            # Adiciona parâmetros da simulação em um dicionário
            parameters[line.split('=')[0]] = (line.split('=')[1]).split('\n')[0]
        elif (line_number + 1) == 6:
            # Pula linha separadora
            pass
        else:
            # Separa valores na lista `x` e na lista `y`
            x.append( float(line.split(',')[0]) )
            y.append( float((line.split(',')[1]).split('\n')[0]) )


def plot_graph(x, y, id):
    methods = {0: 'Upwind', 1:'Superbee', 2:'Van Albada'}
    colors  = {0: 'red',    1: 'green',   2: 'blue'}
    fig,ax  = plt.subplots()
    fig.set_size_inches(12, 7)
    ax.grid(True)
    ax.set(ylim=(0,1))

    plt.suptitle(f"Método {methods[id]}: Concentração X Posição")
    plt.title(rf"$nx = {parameters['nx']}$, " 
              rf"$\Delta t = {parameters['Delta_t']}$, "
              rf"$\Delta x = {parameters['Delta_x']}$, "
              rf"$t_{{final}} = {parameters['t_final']}$, "
              rf"$\bar{{u}} = {parameters['u_bar']}$, ", fontsize=8)

    plt.xlabel("posição x (m)")
    plt.ylabel("concentração Q")
    plt.plot(x,y,'ko-', markerfacecolor=colors[id], markeredgecolor='k')
    plt.draw()

def plot_all_graphs():
    methods = {0: 'Upwind', 1:'Superbee', 2:'Van Albada'}
    colors  = {0: 'red',    1: 'green',   2: 'blue'}
    fig,ax  = plt.subplots()
    fig.set_size_inches(12, 7)
    ax.grid(True)
    ax.set(ylim=(0,1))

    plt.suptitle(f"Todos os métodos: Concentração X Posição")
    plt.title(rf"$nx = {parameters['nx']}$, " 
              rf"$\Delta t = {parameters['Delta_t']}$, "
              rf"$\Delta x = {parameters['Delta_x']}$, "
              rf"$t_{{final}} = {parameters['t_final']}$, "
              rf"$\bar{{u}} = {parameters['u_bar']}$, ", fontsize=8)
    plt.xlabel("posição x (m)")
    plt.ylabel("concentração Q")

    plt.plot(x0, y0, f'{colors[0][0]}o-', markerfacecolor=colors[0], label=methods[0])
    plt.plot(x1, y1, f'{colors[1][0]}o-', markerfacecolor=colors[1], label=methods[1])
    plt.plot(x2, y2, f'{colors[2][0]}o-', markerfacecolor=colors[2], label=methods[2])
    plt.legend(loc='upper right')
    plt.savefig(f"./results/fig/methods_t={parameters['t_final']}.png")
    plt.show()

if __name__ == "__main__":
    main()
