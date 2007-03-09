package network;

import java.io.*;
import java.net.*;
import java.util.*;

import config.*;

/**
 * A <code>Server</code> handles connections from clients, and
 * it sends and recieves data to and from connected clients.
 * Data is handled in such a way the <code>Server</code> doesn't need to know
 * about the content of the data.
 * Whenever a client is connected to a server an <code>InitCommand</code> is
 * send.
 * This command contains data about the scheme that is used in the game and
 * the starting position of this player.
 */
public class Server extends Thread
{
    // The ServerSocket that is used to accept connections.
    private ServerSocket server;

    // Indicates whether this server is running or not.
    private boolean running = true;

    // The list with connected clients
    private LinkedList<Transmitter> clients = new LinkedList<Transmitter>();

    // A list with the recieved Commands from clients that will be
    // send via a GameCommand
    private LinkedList<Command> commands = new LinkedList<Command>();

    // Stream to write log messages to
    private PrintStream log;

    // A unique id that is given to every newly connected client
    // The client_id of the server is Config.SERVER_ID
    private int client_id = 1;

    // The name of the scheme that the players connected to this server have to play.
    private String scheme_name;

    /**
     * Creates a new <code>Server</code> object using a given port number
     * to listen to to accept connections from clients.
     * The <code>Server</code> writes log information to the given
     * <code>PrintStream</code>.
     *
     * @param port        The portnumber for the server to listen on.
     * @param scheme_name The name of the scheme that players connected to
                          this server have to play.
     * @param log         <code>PrintStream</code> to write log messages to.
     * @throws IOException if an I/O error occurs when opening the socket.
     */
    public Server(int port, String scheme_name, PrintStream log)
           throws IOException
    {
        this.scheme_name = scheme_name;
        // Initialize log stream
        this.log = log;
        // Start server
        writeLog("Starting server on port: " + port);
        server = new ServerSocket(port);
        // Start listening for connections
        start();
    }

    /**
     * Listens for connections, and accepts them if possible.
     */
    public void run()
    {
        while (running)
        {
            try
            {
                if (canAccept())
                {
                    // Try to accept a connection from a client.
                    Socket client = server.accept();
                    // Create new transmitter
                    Transmitter n_client = new Transmitter(client, this,
                                                           log, client_id);
                    // Add n_client to connected clients list
                    clients.addLast(n_client);

                    // Send InitCommand to newly connected client.
                    n_client.send(new InitCommand(
                                  Config.SERVER_ID,
                                  client_id,
                                  scheme_name,
                                  (new Random()).nextInt(Config.GAME_WIDTH),
                                  (new Random()).nextInt(Config.GAME_HEIGHT)));

                    // Update client_id, so a client that wants to connects gets
                    // a unique id.
                    client_id++;

                    writeLog("New client accepted!");
                }
                else
                {
                    writeLog("Client refused!");
                }

                // Suspend this Thread for a little while.
                sleep(20);
            }
            catch (InterruptedException ie) {}
            catch (IOException ioe)
            {
                // Thrown when server could
                // not accept connection
                writeLog("Connection failed, message: " + ioe.getMessage() +
                         " (Server closed = " + isClosed() + ").");
            }
        }
    }

    /**
     * Sends a <code>Command</code> to all connected clients.
     *
     * @param command The command to send to all connected clients.
     */
    public void sendToAll(Command command)
    {
        // Send to all clients
        synchronized (clients)
        {
            for (Transmitter t : clients)
            {
                t.send(command);
            }
        }
    }

    /**
     * Recieves a <code>Command</code>, the <code>Server</code> waits
     * for <code>Command</code>s from all other connected clients before
     * they are forwarded in one single <code>GameCommand</code>.
     *
     * @param command The <code>Command</code> command that is recieved.
     */
    public void recieve(Command command)
    {
        writeLog("");
        writeLog("" + command + " recieved from: " + command.getSenderID());
        writeLog("");

        // Add command to command list
        commands.addLast(command);

        // If the commands of all clients are recieved, make a new GameCommand
        if (commands.size() == clients.size())
        {
            // Put commands into an array
            Command[] cs = new Command[commands.size()];
            int i = 0;

            for (Command c : commands)
            {
                cs[i] = c;
                ++i;
            }

            sendToAll(new GameCommand(cs));

            // Clear commands, so a new GameCommand can be sent.
            commands.clear();
        }
    }

    /**
     * Removes a <code>Transmitter</code> from this server.
     * The <code>Transmitter</code> is NOT closed!.
     *
     * @param client The transmitter to remove from the server.
     */
    public void removeClient(Transmitter client)
    {
        clients.remove(client);
        writeLog("Client disconnected.");
    }

    /**
     * Closes all connections to this server, and stops
     * the server from accepting connections.
     */
    public void close()
    {
        running = false;

        // Try to close the server
        try
        {
            server.close();
            writeLog("Server is closed");
        }
        catch (IOException ioe)
        {
            // Thrown if server was blocking on closing
            writeLog("Error while closing server, message " +
                     ioe.getMessage());
        }

        for (Transmitter t : clients)
        {
            t.close();
        }
    }

    /**
     * Returns the closed state of the </code>Server</code>.
     *
     * @Return <code>true</code> if the <code>Server</code> has been closed.
     */
    public boolean isClosed()
    {
        return server.isClosed();
    }

    /**
     * Returns whether it is possible to accept connections from clients.
     *
     * @return <code>true</code> if it is possible to accept connections,
     *         <code>false</code> if this is not possible.
     */
    private boolean canAccept()
    {
        // Only accept when the server is not closed, and there
        // are less than the maximum players connected.
        return !isClosed() &&
               clients.size() < Config.MAX_PLAYERS;
    }

    /**
     * Writes log messages to log if it exsits.
     *
     * @param message The message to write to log.
     */
    private void writeLog(String message)
    {
        if (log != null)
        {
            log.println(message);
        }
    }
}